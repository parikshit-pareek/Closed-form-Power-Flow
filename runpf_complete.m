% MATPOWER CODES COMNBINED TOGETHER IN ONE FILE TO AVOID DEPENDENCIES


function [MVAbase, bus, gen, branch, success, et,J,Ybus] = ...
                runpf_complete(casedata, mpopt, fname, solvedcase)
%RUNPF  Runs a power flow.
%   [RESULTS, SUCCESS] = RUNPF(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs a power flow (full AC Newton's method by default), optionally
%   returning a RESULTS struct and SUCCESS flag.
%
%   Inputs (all are optional):
%       CASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
%       MPOPT : MATPOWER options struct to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more (see also MPOPTION).
%       FNAME : name of a file to which the pretty-printed output will
%           be appended
%       SOLVEDCASE : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       RESULTS : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = succeeded, 0 = failed
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       results = runpf_complete;
%       results = runpf_complete(casedata);
%       results = runpf_complete(casedata, mpopt);
%       results = runpf_complete(casedata, mpopt, fname);
%       results = runpf_complete(casedata, mpopt, fname, solvedcase);
%       [results, success] = runpf_complete(...);
%
%       Alternatively, for compatibility with previous versions of MATPOWER,
%       some of the results can be returned as individual output arguments:
%
%       [baseMVA, bus, gen, branch, success, et] = runpf(...);
%
%   If the pf.enforce_q_lims option is set to true (default is false) then, if
%   any generator reactive power limit is violated after running the AC power
%   flow, the corresponding bus is converted to a PQ bus, with Qg at the
%   limit, and the case is re-run. The voltage magnitude at the bus will
%   deviate from the specified value in order to satisfy the reactive power
%   limit. If the reference bus is converted to PQ, the first remaining PV
%   bus will be used as the slack bus for the next iteration. This may
%   result in the real power output at this generator being slightly off
%   from the specified values.
%
%   Examples:
%       results = runpf('case30');
%       results = runpf('case30', mpoption('pf.enforce_q_lims', 1));
%
%   See also RUNDCPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   Enforcing of generator Q limits inspired by contributions
%   from Mu Lin, Lincoln University, New Zealand (1/14/05).
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default arguments
if nargin < 4
    solvedcase = '';                %% don't save solved case
    if nargin < 3
        fname = '';                 %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casedata = 'case30'; %% default data file is 'case9.m'
            end
        end
    end
end
% mpopt.pf.enforce_q_lims=1;
%% options
qlim = mpopt.pf.enforce_q_lims;         %% enforce Q limits on gens?
 dc = strcmp(upper(mpopt.model), 'DC');  %% use DC formulation?

%% read data
mpc = loadcase(casedata);

%% add zero columns to branch for flows if needed
if size(mpc.branch,2) < QT
  mpc.branch = [ mpc.branch zeros(size(mpc.branch, 1), QT-size(mpc.branch,2)) ];
end

%% convert to internal indexing
mpc = ext2int(mpc);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%%-----  run the power flow  -----
t0 = clock;
its = 0;            %% total iterations
if mpopt.verbose > 0
%     v = mpver('all');
%     fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end
if dc                               %% DC formulation
    if mpopt.verbose > 0
      fprintf(' -- DC Power Flow\n');
    end
    %% initial state
    Va0 = bus(:, VA) * (pi/180);
    
    %% build B matrices and phase shift injections
    [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
    
    %% compute complex bus power injections (generation - load)
    %% adjusted for phase shifters and real shunts
    Pbus = real(makeSbus(baseMVA, bus, gen)) - Pbusinj - bus(:, GS) / baseMVA;
    
    %% "run" the power flow
    [Va, success] = dcpf(B, Pbus, Va0, ref, pv, pq);
    its = 1;
    
    %% update data matrices with solution
    branch(:, [QF, QT]) = zeros(size(branch, 1), 2);
    branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
    branch(:, PT) = -branch(:, PF);
    bus(:, VM) = ones(size(bus, 1), 1);
    bus(:, VA) = Va * (180/pi);
    %% update Pg for slack generator (1st gen at ref bus)
    %% (note: other gens at ref bus are accounted for in Pbus)
    %%      Pg = Pinj + Pload + Gs
    %%      newPg = oldPg + newPinj - oldPinj
    refgen = zeros(size(ref));
    for k = 1:length(ref)
        temp = find(gbus == ref(k));
        refgen(k) = on(temp(1));
    end
    gen(refgen, PG) = gen(refgen, PG) + (B(ref, :) * Va - Pbus(ref)) * baseMVA;
else                                %% AC formulation
    alg = upper(mpopt.pf.alg);
    if mpopt.verbose > 0
        switch alg
            case 'NR'
                solver = 'Newton';
            case 'FDXB'
                solver = 'fast-decoupled, XB';
            case 'FDBX'
                solver = 'fast-decoupled, BX';
            case 'GS'
                solver = 'Gauss-Seidel';
            otherwise
                solver = 'unknown';
        end
%         fprintf(' -- AC Power Flow (%s)\n', solver);
    end
    %% initial state
    % V0    = ones(size(bus, 1), 1);            %% flat start
    V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
    vcb = ones(size(V0));           %% create mask of voltage-controlled buses
    vcb(pq) = 0;                    %% exclude PQ buses
    k = find(vcb(gbus));            %% in-service gens at v-c buses
    V0(gbus(k)) = gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));
    
    if qlim
        ref0 = ref;                         %% save index and angle of
        Varef0 = bus(ref0, VA);             %%   original reference bus(es)
        limited = [];                       %% list of indices of gens @ Q lims
        fixedQg = zeros(size(gen, 1), 1);   %% Qg of gens at Q limits
    end

    %% build admittance matrices
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
    repeat = 1;
    while (repeat)
        %% function for computing V dependent complex bus power injections
        %% (generation - load)
        Sbus = @(Vm)makeSbus(baseMVA, bus, gen, mpopt, Vm);
        
        %% run the power flow
        switch alg
            case 'NR'
                [V, success, iterations,J] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
            case {'FDXB', 'FDBX'}
                [Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
                [V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
            case 'GS'
                if (~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
                        any(mpopt.exp.sys_wide_zip_loads.pw(2:3))) || ...
                        (~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
                        any(mpopt.exp.sys_wide_zip_loads.qw(2:3)))
                    warning('runpf: Gauss-Seidel algorithm does not support ZIP load model. Converting to constant power loads.')
                    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', ...
                                    struct('pw', [], 'qw', []));
                end
                [V, success, iterations] = gausspf(Ybus, Sbus([]), V0, ref, pv, pq, mpopt);
            otherwise
                error('runpf: Only Newton''s method, fast-decoupled, and Gauss-Seidel power flow algorithms currently implemented.');
        end
        its = its + iterations;
        
        %% update data matrices with solution
        [bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
        
        if qlim             %% enforce generator Q limits
            %% find gens with violated Q constraints
            mx = find( gen(:, GEN_STATUS) > 0 ...
                    & gen(:, QG) > gen(:, QMAX) + mpopt.opf.violation );
            mn = find( gen(:, GEN_STATUS) > 0 ...
                    & gen(:, QG) < gen(:, QMIN) - mpopt.opf.violation );
            
            if ~isempty(mx) || ~isempty(mn)  %% we have some Q limit violations
                %% first check for INFEASIBILITY
                infeas = union(mx', mn')';  %% transposes handle fact that
                    %% union of scalars is a row vector
                remaining = find( gen(:, GEN_STATUS) > 0 & ...
                                ( bus(gen(:, GEN_BUS), BUS_TYPE) == PV | ...
                                  bus(gen(:, GEN_BUS), BUS_TYPE) == REF ));
                if length(infeas) == length(remaining) && all(infeas == remaining) && ...
                        (isempty(mx) || isempty(mn))
                    %% all remaining PV/REF gens are violating AND all are
                    %% violating same limit (all violating Qmin or all Qmax)
                    if mpopt.verbose
                        fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(infeas));
                    end
                    success = 0;
                    break;
                end

                %% one at a time?
                if qlim == 2    %% fix largest violation, ignore the rest
                    [junk, k] = max([gen(mx, QG) - gen(mx, QMAX);
                                     gen(mn, QMIN) - gen(mn, QG)]);
                    if k > length(mx)
                        mn = mn(k-length(mx));
                        mx = [];
                    else
                        mx = mx(k);
                        mn = [];
                    end
                end

                if mpopt.verbose && ~isempty(mx)
                    fprintf('Gen %d at upper Q limit, converting to PQ bus\n', mx);
                end
                if mpopt.verbose && ~isempty(mn)
                    fprintf('Gen %d at lower Q limit, converting to PQ bus\n', mn);
                end
                
                %% save corresponding limit values
                fixedQg(mx) = gen(mx, QMAX);
                fixedQg(mn) = gen(mn, QMIN);
                mx = [mx;mn];
                
                %% convert to PQ bus
                gen(mx, QG) = fixedQg(mx);      %% set Qg to binding limit
                gen(mx, GEN_STATUS) = 0;        %% temporarily turn off gen,
                for i = 1:length(mx)            %% (one at a time, since
                    bi = gen(mx(i), GEN_BUS);   %%  they may be at same bus)
                    bus(bi, [PD,QD]) = ...      %% adjust load accordingly,
                        bus(bi, [PD,QD]) - gen(mx(i), [PG,QG]);
                end
                if length(ref) > 1 && any(bus(gen(mx, GEN_BUS), BUS_TYPE) == REF)
                    error('runpf: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                end
                bus(gen(mx, GEN_BUS), BUS_TYPE) = PQ;   %% & set bus type to PQ
                
                %% update bus index lists of each type of bus
                ref_temp = ref;
                [ref, pv, pq] = bustypes(bus, gen);
                %% previous line can modify lists to select new REF bus
                %% if there was none, so we should update bus with these
                %% just to keep them consistent
                if ref ~= ref_temp
                    bus(ref, BUS_TYPE) = REF;
                    bus( pv, BUS_TYPE) = PV;
                    if mpopt.verbose
                        fprintf('Bus %d is new slack bus\n', ref);
                    end
                end
                limited = [limited; mx];
            else
                repeat = 0; %% no more generator Q limits violated
            end
        else
            repeat = 0;     %% don't enforce generator Q limits, once is enough
        end
    end
    if qlim && ~isempty(limited)
        %% restore injections from limited gens (those at Q limits)
        gen(limited, QG) = fixedQg(limited);    %% restore Qg value,
        for i = 1:length(limited)               %% (one at a time, since
            bi = gen(limited(i), GEN_BUS);      %%  they may be at same bus)
            bus(bi, [PD,QD]) = ...              %% re-adjust load,
                bus(bi, [PD,QD]) + gen(limited(i), [PG,QG]);
        end
        gen(limited, GEN_STATUS) = 1;               %% and turn gen back on
        if ref ~= ref0
            %% adjust voltage angles to make original ref bus correct
            bus(:, VA) = bus(:, VA) - bus(ref0, VA) + Varef0;
        end
    end
end
mpc.et = etime(clock, t0);
mpc.success = success;
mpc.iterations = its;

%%-----  output results  -----
%% convert back to original bus numbering & print results
[mpc.bus, mpc.gen, mpc.branch] = deal(bus, gen, branch);
results = int2ext(mpc);

%% zero out result fields of out-of-service gens & branches
if ~isempty(results.order.gen.status.off)
  results.gen(results.order.gen.status.off, [PG QG]) = 0;
end
if ~isempty(results.order.branch.status.off)
  results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
end

if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        if mpopt.out.all == 0
%             printpf(results, fd, mpoption(mpopt, 'out.all', -1));
        else
%             printpf(results, fd, mpopt);
        end
        fclose(fd);
    end
end

% printpf(results, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, results);
end

if nargout == 1 || nargout == 2
    MVAbase = results;
    bus = success;
elseif nargout > 2
    [MVAbase, bus, gen, branch, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, results.et);
% else  %% don't define MVAbase, so it doesn't print anything
end


%% Function involved in this 

function varargout = deal(varargin)
%DEAL Deal inputs to outputs.
%    [A,B,C,...] = DEAL(X,Y,Z,...) simply matches up the input and
%       output lists.  It is the same as A=X, B=Y, C=Z, ...
%    [A,B,C,...] = DEAL(X) copies the single input to all
%       the requested outputs.  It is the same as A=X, B=X, C=X, ...
%
%    DEAL is most useful when used with cell arrays and structures
%    via comma separated list expansion.  Here are some useful
%    constructions:
%    [S.FIELD] = DEAL(X) sets all the fields with the name FIELD
%       in the structure array S to the value X.  If S doesn't
%       exist, use [S(1:M).FIELD] = DEAL(X);
%    [X{:}] = DEAL(A.FIELD) copies the values of the field with
%       name FIELD to the cell array X.  If X doesn't exist,
%       use [X{1:M}] = DEAL(A.FIELD).
%    [A,B,C,...] = DEAL(X{:}) copies the contents of the cell
%       array X to the separate variables A,B,C,...
%    [A,B,C,...] = DEAL(S.FIELD) copies the contents of the fields
%       with the name FIELD to separate variables A,B,C,...
%
%    Examples:
%       sys = {rand(3) ones(3,1) eye(3) zeros(3,1)};
%       [a,b,c,d] = deal(sys{:});
%
%       direc = dir; filenames = {};
%       [filenames{1:length(direc),1}] = deal(direc.name);
%
%    See also LISTS, PAREN.

%   Copyright 1984-2005 The MathWorks, Inc. 

if nargin==1,
  varargout = varargin(ones(1,nargout));
else
  if nargout ~= nargin
    error(message('MATLAB:deal:narginNargoutMismatch'))
  end
  varargout = varargin;
end

function [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus
%IDX_BUS   Defines constants for named column indices to bus matrix.
%   Example:
%
%   [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%   VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pd = bus(4, PD);       % get the real power demand at bus 4
%    bus(:, VMIN) = 0.95;   % set the min voltage magnitude to 0.95 at all buses
% 
%   The index, name and meaning of each column of the bus matrix is given
%   below:
%
%   columns 1-13 must be included in input matrix (in case file)
%    1  BUS_I       bus number (positive integer)
%    2  BUS_TYPE    bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
%    3  PD          Pd, real power demand (MW)
%    4  QD          Qd, reactive power demand (MVAr)
%    5  GS          Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
%    6  BS          Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
%    7  BUS_AREA    area number, (positive integer)
%    8  VM          Vm, voltage magnitude (p.u.)
%    9  VA          Va, voltage angle (degrees)
%    10 BASE_KV     baseKV, base voltage (kV)
%    11 ZONE        zone, loss zone (positive integer)
%    12 VMAX        maxVm, maximum voltage magnitude (p.u.)
%    13 VMIN        minVm, minimum voltage magnitude (p.u.)
%   
%   columns 14-17 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    14 LAM_P       Lagrange multiplier on real power mismatch (u/MW)
%    15 LAM_Q       Lagrange multiplier on reactive power mismatch (u/MVAr)
%    16 MU_VMAX     Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
%    17 MU_VMIN     Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
% 
%   additional constants, used to assign/compare values in the BUS_TYPE column
%    1  PQ    PQ bus
%    2  PV    PV bus
%    3  REF   reference bus
%    4  NONE  isolated bus
%
%   See also DEFINE_CONSTANTS.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define bus types
PQ      = 1;
PV      = 2;
REF     = 3;
NONE    = 4;

%% define the indices
BUS_I       = 1;    %% bus number (1 to 29997)
BUS_TYPE    = 2;    %% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD          = 3;    %% Pd, real power demand (MW)
QD          = 4;    %% Qd, reactive power demand (MVAr)
GS          = 5;    %% Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 6;    %% Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
BUS_AREA    = 7;    %% area number, 1-100
VM          = 8;    %% Vm, voltage magnitude (p.u.)
VA          = 9;    %% Va, voltage angle (degrees)
BASE_KV     = 10;   %% baseKV, base voltage (kV)
ZONE        = 11;   %% zone, loss zone (1-999)
VMAX        = 12;   %% maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
VMIN        = 13;   %% minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
LAM_P       = 14;   %% Lagrange multiplier on real power mismatch (u/MW)
LAM_Q       = 15;   %% Lagrange multiplier on reactive power mismatch (u/MVAr)
MU_VMAX     = 16;   %% Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN     = 17;   %% Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)


function [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen
%IDX_GEN   Defines constants for named column indices to gen matrix.
%   Example:
%
%   [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
%   MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
%   QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pg = gen(4, PG);   % get the real power output of generator 4
%    gen(:, PMIN) = 0;  % set to zero the minimum real power limit of all gens
% 
%   The index, name and meaning of each column of the gen matrix is given
%   below:
%
%   columns 1-21 must be included in input matrix (in case file)
%    1  GEN_BUS     bus number
%    2  PG          Pg, real power output (MW)
%    3  QG          Qg, reactive power output (MVAr)
%    4  QMAX        Qmax, maximum reactive power output (MVAr)
%    5  QMIN        Qmin, minimum reactive power output (MVAr)
%    6  VG          Vg, voltage magnitude setpoint (p.u.)
%    7  MBASE       mBase, total MVA base of machine, defaults to baseMVA
%    8  GEN_STATUS  status, > 0 - in service, <= 0 - out of service
%    9  PMAX        Pmax, maximum real power output (MW)
%    10 PMIN        Pmin, minimum real power output (MW)
%    11 PC1         Pc1, lower real power output of PQ capability curve (MW)
%    12 PC2         Pc2, upper real power output of PQ capability curve (MW)
%    13 QC1MIN      Qc1min, minimum reactive power output at Pc1 (MVAr)
%    14 QC1MAX      Qc1max, maximum reactive power output at Pc1 (MVAr)
%    15 QC2MIN      Qc2min, minimum reactive power output at Pc2 (MVAr)
%    16 QC2MAX      Qc2max, maximum reactive power output at Pc2 (MVAr)
%    17 RAMP_AGC    ramp rate for load following/AGC (MW/min)
%    18 RAMP_10     ramp rate for 10 minute reserves (MW)
%    19 RAMP_30     ramp rate for 30 minute reserves (MW)
%    20 RAMP_Q      ramp rate for reactive power (2 sec timescale) (MVAr/min)
%    21 APF         area participation factor
%   
%   columns 22-25 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    22 MU_PMAX     Kuhn-Tucker multiplier on upper Pg limit (u/MW)
%    23 MU_PMIN     Kuhn-Tucker multiplier on lower Pg limit (u/MW)
%    24 MU_QMAX     Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
%    25 MU_QMIN     Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)
%
%   See also DEFINE_CONSTANTS.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
GEN_BUS     = 1;    %% bus number
PG          = 2;    %% Pg, real power output (MW)
QG          = 3;    %% Qg, reactive power output (MVAr)
QMAX        = 4;    %% Qmax, maximum reactive power output at Pmin (MVAr)
QMIN        = 5;    %% Qmin, minimum reactive power output at Pmin (MVAr)
VG          = 6;    %% Vg, voltage magnitude setpoint (p.u.)
MBASE       = 7;    %% mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 8;    %% status, 1 - machine in service, 0 - machine out of service
PMAX        = 9;    %% Pmax, maximum real power output (MW)
PMIN        = 10;   %% Pmin, minimum real power output (MW)
PC1         = 11;   %% Pc1, lower real power output of PQ capability curve (MW)
PC2         = 12;   %% Pc2, upper real power output of PQ capability curve (MW)
QC1MIN      = 13;   %% Qc1min, minimum reactive power output at Pc1 (MVAr)
QC1MAX      = 14;   %% Qc1max, maximum reactive power output at Pc1 (MVAr)
QC2MIN      = 15;   %% Qc2min, minimum reactive power output at Pc2 (MVAr)
QC2MAX      = 16;   %% Qc2max, maximum reactive power output at Pc2 (MVAr)
RAMP_AGC    = 17;   %% ramp rate for load following/AGC (MW/min)
RAMP_10     = 18;   %% ramp rate for 10 minute reserves (MW)
RAMP_30     = 19;   %% ramp rate for 30 minute reserves (MW)
RAMP_Q      = 20;   %% ramp rate for reactive power (2 sec timescale) (MVAr/min)
APF         = 21;   %% area participation factor

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_PMAX     = 22;   %% Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN     = 23;   %% Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX     = 24;   %% Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN     = 25;   %% Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

%% Note: When a generator's PQ capability curve is not simply a box and the
%% upper Qg limit is binding, the multiplier on this constraint is split into
%% it's P and Q components and combined with the appropriate MU_Pxxx and
%% MU_Qxxx values. Likewise for the lower Q limits.

function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch
%IDX_BRCH   Defines constants for named column indices to branch matrix.
%   Example:
%
%   [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%   TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%   ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    branch(4, BR_STATUS) = 0;              % take branch 4 out of service
%    Ploss = branch(:, PF) + branch(:, PT); % compute real power loss vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-11 must be included in input matrix (in case file)
%    1  F_BUS       f, from bus number
%    2  T_BUS       t, to bus number
%    3  BR_R        r, resistance (p.u.)
%    4  BR_X        x, reactance (p.u.)
%    5  BR_B        b, total line charging susceptance (p.u.)
%    6  RATE_A      rateA, MVA rating A (long term rating)
%    7  RATE_B      rateB, MVA rating B (short term rating)
%    8  RATE_C      rateC, MVA rating C (emergency rating)
%    9  TAP         ratio, transformer off nominal turns ratio
%    10 SHIFT       angle, transformer phase shift angle (degrees)
%    11 BR_STATUS   initial branch status, 1 - in service, 0 - out of service
%    12 ANGMIN      minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%    13 ANGMAX      maximum angle difference, angle(Vf) - angle(Vt) (degrees)
%                   (The voltage angle difference is taken to be unbounded below
%                    if ANGMIN < -360 and unbounded above if ANGMAX > 360.
%                    If both parameters are zero, it is unconstrained.)
%
%   columns 14-17 are added to matrix after power flow or OPF solution
%   they are typically not present in the input matrix
%    14 PF          real power injected at "from" bus end (MW)
%    15 QF          reactive power injected at "from" bus end (MVAr)
%    16 PT          real power injected at "to" bus end (MW)
%    17 QT          reactive power injected at "to" bus end (MVAr)
%
%   columns 18-21 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    18 MU_SF       Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
%    19 MU_ST       Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
%    20 MU_ANGMIN   Kuhn-Tucker multiplier lower angle difference limit (u/degree)
%    21 MU_ANGMAX   Kuhn-Tucker multiplier upper angle difference limit (u/degree)
%
%   See also DEFINE_CONSTANTS.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
F_BUS       = 1;    %% f, from bus number
T_BUS       = 2;    %% t, to bus number
BR_R        = 3;    %% r, resistance (p.u.)
BR_X        = 4;    %% x, reactance (p.u.)
BR_B        = 5;    %% b, total line charging susceptance (p.u.)
RATE_A      = 6;    %% rateA, MVA rating A (long term rating)
RATE_B      = 7;    %% rateB, MVA rating B (short term rating)
RATE_C      = 8;    %% rateC, MVA rating C (emergency rating)
TAP         = 9;    %% ratio, transformer off nominal turns ratio
SHIFT       = 10;   %% angle, transformer phase shift angle (degrees)
BR_STATUS   = 11;   %% initial branch status, 1 - in service, 0 - out of service
ANGMIN      = 12;   %% minimum angle difference, angle(Vf) - angle(Vt) (degrees)
ANGMAX      = 13;   %% maximum angle difference, angle(Vf) - angle(Vt) (degrees)

%% included in power flow solution, not necessarily in input
PF          = 14;   %% real power injected at "from" bus end (MW)       (not in PTI format)
QF          = 15;   %% reactive power injected at "from" bus end (MVAr) (not in PTI format)
PT          = 16;   %% real power injected at "to" bus end (MW)         (not in PTI format)
QT          = 17;   %% reactive power injected at "to" bus end (MVAr)   (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_SF       = 18;   %% Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 19;   %% Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
MU_ANGMIN   = 20;   %% Kuhn-Tucker multiplier lower angle difference limit (u/degree)
MU_ANGMAX   = 21;   %% Kuhn-Tucker multiplier upper angle difference limit (u/degree)













function [ref, pv, pq] = bustypes(bus, gen)
%BUSTYPES   Builds index lists for each type of bus (REF, PV, PQ).
%   [REF, PV, PQ] = BUSTYPES(BUS, GEN)
%   Generators with "out-of-service" status are treated as PQ buses with
%   zero generation (regardless of Pg/Qg values in gen). Expects BUS and
%   GEN have been converted to use internal consecutive bus numbering.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get generator status
% bus_gen_status = zeros(size(bus, 1), 1);
% bus_gen_status(gen(:, GEN_BUS)) = gen(:, GEN_STATUS) > 0;
nb = size(bus, 1);
ng = size(gen, 1);
Cg = sparse(gen(:, GEN_BUS), (1:ng)', gen(:, GEN_STATUS) > 0, nb, ng);  %% gen connection matrix
                                        %% element i, j is 1 if, generator j at bus i is ON
bus_gen_status = Cg * ones(ng, 1);      %% number of generators at each bus that are ON


%% form index lists for slack, PV, and PQ buses
ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   %% reference bus index
pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   %% PV bus indices
pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   %% PQ bus indices

%% pick a new reference bus if for some reason there is none (may have been shut down)
if isempty(ref)
    if isempty(pv)
        %% no PV bus left to convert to reference bus
    else
        ref = pv(1);    %% use the first PV bus
        pv(1) = [];     %% delete it from PV list
    end
end

function [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch)
%MAKEYBUS   Builds the bus admittance matrix and branch admittance matrices.
%   [YBUS, YF, YT] = MAKEYBUS(MPC)
%   [YBUS, YF, YT] = MAKEYBUS(BASEMVA, BUS, BRANCH)
%   
%   Returns the full bus admittance matrix (i.e. for all buses) and the
%   matrices YF and YT which, when multiplied by a complex voltage vector,
%   yield the vector currents injected into each line from the "from" and
%   "to" buses respectively of each line. Does appropriate conversions to p.u.
%   Inputs can be a MATPOWER case struct or individual BASEMVA, BUS and
%   BRANCH values. Bus numbers must be consecutive beginning at 1
%   (i.e. internal ordering).
%
%   See also MAKEJAC, MAKESBUS, EXT2INT.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeYbus: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%% for each branch, compute the elements of the branch admittance matrix where
%%
%%      | If |   | Yff  Yft |   | Vf |
%%      |    | = |          | * |    |
%%      | It |   | Ytf  Ytt |   | Vt |
%%
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
Bc = stat .* branch(:, BR_B);                           %% line charging susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% compute shunt admittance
%% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
%% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
%% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
%% i.e. Ysh = Psh + j Qsh, so ...
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA; %% vector of shunt admittances

%% build connection matrices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses

%% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
%% at each branch's "from" bus, and Yt is the same for the "to" bus end
i = [1:nl; 1:nl]';                              %% double set of row indices
Yf = sparse(i, [f; t], [Yff; Yft], nl, nb);
Yt = sparse(i, [f; t], [Ytf; Ytt], nl, nb);
% Yf = spdiags(Yff, 0, nl, nl) * Cf + spdiags(Yft, 0, nl, nl) * Ct;
% Yt = spdiags(Ytf, 0, nl, nl) * Cf + spdiags(Ytt, 0, nl, nl) * Ct;

%% build Ybus
Ybus = Cf' * Yf + Ct' * Yt + ...                %% branch admittances
        sparse(1:nb, 1:nb, Ysh, nb, nb);        %% shunt admittance

    function [Sbus, dSbus_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm, Sg)
%MAKESBUS   Builds the vector of complex bus power injections.
%   SBUS = MAKESBUS(BASEMVA, BUS, GEN)
%   SBUS = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM)
%   SBUS = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM, SG)
%   returns the vector of complex bus power injections, that is, generation
%   minus load. Power is expressed in per unit. If the MPOPT and VM arguments
%   are present it evaluates any ZIP loads based on the provided voltage
%   magnitude vector. If VM is empty, it assumes nominal voltage. If SG is
%   provided, it is a complex ng x 1 vector of generator power injections in
%   p.u., and overrides the PG and QG columns in GEN, using GEN only for
%   connectivity information.
%
%   [SBUS, DSBUS_DVM] = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM)
%   With two output arguments, it computes the partial derivative of the
%   bus injections with respect to voltage magnitude, leaving the first
%   return value SBUS empty. If VM is empty, it assumes no voltage dependence
%   and returns a sparse zero matrix.
%
%   See also MAKEYBUS.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default inputs
if nargin < 5
    Vm = [];
    if nargin < 4
        mpopt = [];
    end
end
nb = size(bus, 1);

%% get load parameters
Sd = makeSdzip(baseMVA, bus, mpopt);

if nargout == 2
    Sbus = [];
    if isempty(Vm)
        dSbus_dVm = sparse(nb, nb);
    else
        dSbus_dVm = -(spdiags(Sd.i + 2 * Vm .* Sd.z, 0, nb, nb));
    end
else
    %% compute per-bus generation in p.u.
    on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
    gbus = gen(on, GEN_BUS);                %% what buses are they at?
    ngon = size(on, 1);
    Cg = sparse(gbus, (1:ngon)', 1, nb, ngon);  %% connection matrix
                                                %% element i, j is 1 if
                                                %% gen on(j) at bus i is ON
    if nargin > 5 && ~isempty(Sg)
        Sbusg = Cg * Sg(on);
    else
        Sbusg = Cg * (gen(on, PG) + 1j * gen(on, QG)) / baseMVA;
    end

    %% compute per-bus loads in p.u.
    if isempty(Vm)
        Vm = ones(nb, 1);
    end
    Sbusd = Sd.p + Sd.i .* Vm + Sd.z .* Vm.^2;

    %% form net complex bus power injection vector
    %% (power injected by generators + power injected by loads)
    Sbus = Sbusg - Sbusd;
end

function [V, converged, i,J] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves the power flow using a full Newton's method.
%   [V, CONVERGED, I] = NEWTONPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. MPOPT is a MATPOWER options struct which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see MPOPTION for details). Uses default options if
%   this parameter is not given. Returns the final complex voltages, a
%   flag which indicates whether it converged or not, and the number of
%   iterations performed.
%
%   See also RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
mpopt.pf.nr.max_it=20;
%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.nr.max_it;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus(Vm);
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% check tolerance
normF = norm(F, inf);
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% evaluate Jacobian
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
    [dummy, neg_dSd_dVm] = Sbus(Vm);
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
    
    J = [   j11 j12;
            j21 j22;    ];

    %% compute update step
    dx = -(J \ F);

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    F = [   real(mis(pv));
            real(mis(pq));
            imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
%             fprintf('\nNewton''s method power flow converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow did not converge in %d iterations.\n', i);
    end
end

function [Bp, Bpp] = makeB(baseMVA, bus, branch, alg)
%MAKEB   Builds the FDPF matrices, B prime and B double prime.
%   [BP, BPP] = MAKEB(MPC, ALG)
%   [BP, BPP] = MAKEB(BASEMVA, BUS, BRANCH, ALG)
%
%   Returns the two matrices B prime and B double prime used in the fast
%   decoupled power flow. Does appropriate conversions to p.u. ALG is either
%   'FDXB' or 'FDBX', the corresponding value of MPOPT.pf.alg option
%   specifying the power flow algorithm.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   Note: For backward compatibility, ALG can also take on a value of
%   2 or 3, corresponding to values of the old PF_ALG option. This usage
%   is deprecated and will be removed in a future version.
%
%   Example:
%       [Bp, Bpp] = makeB(baseMVA, bus, branch, 'FDXB');
%
%   See also FDPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    if nargin == 2
        alg = bus;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% backward compatiblility (deprecated)
if ~ischar(alg)
    if alg == 2
        alg = 'FDXB';
    elseif alg == 3
        alg = 'FDBX';
    end
end

%% check for valid ALG value
alg = upper(alg);
if ~strcmp(alg, 'FDXB') && ~strcmp(alg, 'FDBX')
    error('makeB: ''%s'' is not a valid value for ALG', alg);
end

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_bus = bus;                             %% modify a copy of bus
temp_bus(:, BS) = zeros(nb, 1);             %% zero out shunts at buses
temp_branch(:, BR_B) = zeros(nl, 1);        %% zero out line charging shunts
temp_branch(:, TAP) = ones(nl, 1);          %% cancel out taps
if strcmp(alg, 'FDXB')                      %% if XB method
    temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
end
Bp = -imag( makeYbus(baseMVA, temp_bus, temp_branch) );

%%-----  form Bpp (B double prime)  -----
if nargout == 2
    temp_branch = branch;                       %% modify a copy of branch
    temp_branch(:, SHIFT) = zeros(nl, 1);       %% zero out phase shifters
    if strcmp(alg, 'FDBX')                      %% if BX method
        temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
    end
    Bpp = -imag( makeYbus(baseMVA, bus, temp_branch) );
end

function [V, converged, i] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%FDPF  Solves the power flow using a fast decoupled method.
%   [V, CONVERGED, I] = FDPF(YBUS, SBUS, V0, BP, BPP, REF, PV, PQ, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, the FDPF matrices B prime
%   and B double prime, and column vectors with the lists of bus indices
%   for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
%   vector contains the set point for generator (including ref bus)
%   buses, and the reference angle of the swing bus, as well as an initial
%   guess for remaining magnitudes and angles. MPOPT is a MATPOWER options
%   vector which can be used to set the termination tolerance, maximum
%   number of iterations, and output options (see MPOPTION for details).
%   Uses default options if this parameter is not given. Returns the
%   final complex voltages, a flag which indicates whether it converged
%   or not, and the number of iterations performed.
%
%   See also RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.fd.max_it;
if have_fcn('matlab') && have_fcn('matlab', 'vnum') < 7.3
    lu_vec = 0;     %% lu(..., 'vector') syntax not supported
else
    lu_vec = 1;
end

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);

%% evaluate initial mismatch
mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
P = real(mis([pv; pq]));
Q = imag(mis(pq));

%% check tolerance
normP = norm(P, inf);
normQ = norm(Q, inf);
if mpopt.verbose > 1
    fprintf('\niteration     max mismatch (p.u.)  ');
    fprintf('\ntype   #        P            Q     ');
    fprintf('\n---- ----  -----------  -----------');
    fprintf('\n  -  %3d   %10.3e   %10.3e', i, normP, normQ);
end
if normP < tol && normQ < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% reduce B matrices
Bp = Bp([pv; pq], [pv; pq]);
Bpp = Bpp(pq, pq);

%% factor B matrices
if lu_vec
    [Lp,  Up,  pp,  qp ] = lu(Bp,  'vector');
    [Lpp, Upp, ppp, qpp] = lu(Bpp, 'vector');
    [junk, iqp ] = sort(qp);
    [junk, iqpp] = sort(qpp);
    % [~, iqp ] = sort(qp);
    % [~, iqpp] = sort(qpp);
else
    [Lp, Up, Pp] = lu(Bp);
    [Lpp, Upp, Ppp] = lu(Bpp);
end

%% do P and Q iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %%-----  do P iteration, update Va  -----
    if lu_vec
        dVa = -( Up \  (Lp \ P(pp)) );
        dVa = dVa(iqp);
    else
        dVa = -( Up \  (Lp \ (Pp * P)));
    end

    %% update voltage
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);
    if mpopt.verbose > 1
        fprintf('\n  P  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i-1);
        end
        break;
    end

    %%-----  do Q iteration, update Vm  -----
    if lu_vec
        dVm = -( Upp \ (Lpp \ Q(ppp)) );
        dVm = dVm(iqpp);
    else
        dVm = -( Upp \ (Lpp \ (Ppp * Q)) );
    end

    %% update voltage
    Vm(pq) = Vm(pq) + dVm;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);
    if mpopt.verbose > 1
        fprintf('\n  Q  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i);
        end
        break;
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nFast-decoupled power flow did not converge in %d iterations.\n', i);
    end
end

function [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq, mpopt)
%PFSOLN  Updates bus, gen, branch data structures to match power flow soln.
%   [BUS, GEN, BRANCH] = PFSOLN(BASEMVA, BUS0, GEN0, BRANCH0, ...
%                                   YBUS, YF, YT, V, REF, PV, PQ, MPOPT)

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default options
if nargin < 12
    mpopt = mpoption();
end

%% initialize return values
bus     = bus0;
gen     = gen0;
branch  = branch0;

%%----- update bus voltages -----
bus(:, VM) = abs(V);
bus(:, VA) = angle(V) * 180 / pi;

%%----- update Qg for gens at PV/slack buses and Pg for slack bus(es) -----
%% generator info
on = find(gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
        bus(gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
off = find(gen(:, GEN_STATUS) <= 0);    %% which generators are off?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%% compute total injected bus powers
Sbus = V(gbus) .* conj(Ybus(gbus, :) * V);

%% update Qg for generators at PV/slack buses
gen(off, QG) = zeros(length(off), 1);   %% zero out off-line Qg
%% don't touch the ones at PQ buses
[Pd_gbus, Qd_gbus] = total_load(bus(gbus, :), [], 'bus', [], mpopt);
gen(on, QG) = imag(Sbus) * baseMVA + Qd_gbus;   %% inj Q + local Qd
%% ... at this point any buses with more than one generator will have
%% the total Q dispatch for the bus assigned to each generator. This
%% must be split between them. We do it first equally, then in proportion
%% to the reactive range of the generator.

if length(on) > 1
    %% build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
    nb = size(bus, 1);
    ngon = size(on, 1);
    Cg = sparse((1:ngon)', gbus, ones(ngon, 1), ngon, nb);

    %% divide Qg by number of generators at the bus to distribute equally
    ngg = Cg * sum(Cg)';    %% ngon x 1, number of gens at this gen's bus
    gen(on, QG) = gen(on, QG) ./ ngg;

    %% divide proportionally
    Cmin = sparse((1:ngon)', gbus, gen(on, QMIN), ngon, nb);
    Cmax = sparse((1:ngon)', gbus, gen(on, QMAX), ngon, nb);
    Qg_tot = Cg' * gen(on, QG);     %% nb x 1 vector of total Qg at each bus
    Qg_min = sum(Cmin)';            %% nb x 1 vector of min total Qg at each bus
    Qg_max = sum(Cmax)';            %% nb x 1 vector of max total Qg at each bus
    ig = find(Cg * Qg_min == Cg * Qg_max);  %% gens at buses with Qg range = 0
    Qg_save = gen(on(ig), QG);
    gen(on, QG) = gen(on, QMIN) + ...
        (Cg * ((Qg_tot - Qg_min)./(Qg_max - Qg_min + eps))) .* ...
            (gen(on, QMAX) - gen(on, QMIN));    %%    ^ avoid div by 0
    gen(on(ig), QG) = Qg_save;
end                                             %% (terms are mult by 0 anyway)

%% update Pg for slack gen(s)
for k = 1:length(ref)
    refgen = find(gbus == ref(k));              %% which is(are) the reference gen(s)?
    Pd_refk = total_load(bus(ref(k), :), [], 'bus', [], mpopt);
    gen(on(refgen(1)), PG) = real(Sbus(refgen(1))) * baseMVA + Pd_refk; %% inj P + local Pd
    if length(refgen) > 1       %% more than one generator at this ref bus
        %% subtract off what is generated by other gens at this bus
        gen(on(refgen(1)), PG) = gen(on(refgen(1)), PG) ...
                                - sum(gen(on(refgen(2:length(refgen))), PG));
    end
end

%%----- update/compute branch power flows -----
out = find(branch(:, BR_STATUS) == 0);      %% out-of-service branches
br = find(branch(:, BR_STATUS));            %% in-service branches
Sf = V(branch(br, F_BUS)) .* conj(Yf(br, :) * V) * baseMVA; %% complex power at "from" bus
St = V(branch(br, T_BUS)) .* conj(Yt(br, :) * V) * baseMVA; %% complex power injected at "to" bus
branch(br, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];
branch(out, [PF, QF, PT, QT]) = zeros(length(out), 4);



function opt = mpoption(varargin)
%MPOPTION  Used to set and retrieve a MATPOWER options struct.
%
%   OPT = MPOPTION
%       Returns the default options struct.
%
%   OPT = MPOPTION(OVERRIDES)
%       Returns the default options struct, with some fields overridden
%       by values from OVERRIDES, which can be a struct or the name of
%       a function that returns a struct.
%
%   OPT = MPOPTION(NAME1, VALUE1, NAME2, VALUE2, ...)
%       Same as previous, except override options are specified by NAME,
%       VALUE pairs. This can be used to set any part of the options
%       struct. The names can be individual fields or multi-level field
%       names with embedded periods. The values can be scalars or structs.
%
%       For backward compatibility, the NAMES and VALUES may correspond
%       to old-style MATPOWER option names (elements in the old-style
%       options vector) as well.
%
%   OPT = MPOPTION(OPT0)
%       Converts an old-style options vector OPT0 into the corresponding
%       options struct. If OPT0 is an options struct it does nothing.
%
%   OPT = MPOPTION(OPT0, OVERRIDES)
%       Applies overrides to an existing set of options, OPT0, which
%       can be an old-style options vector or an options struct.
%
%   OPT = MPOPTION(OPT0, NAME1, VALUE1, NAME2, VALUE2, ...)
%       Same as above except it uses the old-style options vector OPT0
%       as a base instead of the old default options vector.
%
%   OPT_VECTOR = MPOPTION(OPT, [])
%       Creates and returns an old-style options vector from an
%       options struct OPT.
%
%   Note: The use of old-style MATPOWER options vectors and their
%         names and values has been deprecated and will be removed
%         in a future version of MATPOWER. Until then, all uppercase
%         option names are not permitted for new top-level options.
%
%   Examples:
%       mpopt = mpoption('pf.alg', 'FDXB', 'pf.tol', 1e-4);
%       mpopt = mpoption(mpopt, 'opf.dc.solver', 'CPLEX', 'verbose', 2);
%
%The currently defined options are as follows:
%
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%Model options:
%   model                   'AC'        AC vs. DC power flow model
%       [ 'AC' - use nonlinear AC model & corresponding algorithms/options  ]
%       [ 'DC' - use linear DC model & corresponding algorithms/options     ]
%
%Power Flow options:
%   pf.alg                  'NR'        AC power flow algorithm
%       [ 'NR'   - Newton's method                                          ]
%       [ 'FDXB' - Fast-Decoupled (XB version)                              ]
%       [ 'FDBX' - Fast-Decoupled (BX version)                              ]
%       [ 'GS'   - Gauss-Seidel                                             ]
%   pf.tol                  1e-8        termination tolerance on per unit
%                                       P & Q mismatch
%   pf.nr.max_it            10          maximum number of iterations for
%                                       Newton's method
%   pf.fd.max_it            30          maximum number of iterations for
%                                       fast decoupled method
%   pf.gs.max_it            1000        maximum number of iterations for
%                                       Gauss-Seidel method
%   pf.enforce_q_lims       0           enforce gen reactive power limits at
%                                       expense of |V|
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%       [  2 - enforce limits, one-at-a-time bus type conversion            ]
%
%Continuation Power Flow options:
%   cpf.parameterization    3           choice of parameterization
%       [  1 - natural                                                      ]
%       [  2 - arc length                                                   ]
%       [  3 - pseudo arc length                                            ]
%   cpf.stop_at             'NOSE'      determins stopping criterion
%       [ 'NOSE'     - stop when nose point is reached                      ]
%       [ 'FULL'     - trace full nose curve                                ]
%       [ <lam_stop> - stop upon reaching specified target lambda value     ]
%   cpf.enforce_p_lims      0           enforce gen active power limits
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%   cpf.enforce_q_lims      0           enforce gen reactive power limits at
%                                       expense of |V|
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%   cpf.step                0.05        continuation power flow step size
%   cpf.adapt_step          0           toggle adaptive step size feature
%       [  0 - adaptive step size disabled                                  ]
%       [  1 - adaptive step size enabled                                   ]
%   cpf.step_min            1e-4        minimum allowed step size
%   cpf.step_max            0.2         maximum allowed step size
%   cpf.adapt_step_damping  0.7         damping factor for adaptive step
%                                       sizing
%   cpf.adapt_step_tol      1e-3        tolerance for adaptive step sizing
%   cpf.target_lam_tol      1e-5        tolerance for target lambda detection
%   cpf.nose_tol            1e-5        tolerance for nose point detection (pu)
%   cpf.p_lims_tol          0.01        tolerance for generator active
%                                       power limit enforcement (MW)
%   cpf.q_lims_tol          0.01        tolerance for generator reactive
%                                       power limit enforcement (MVAR)
%   cpf.plot.level          0           control plotting of noze curve
%       [  0 - do not plot nose curve                                       ]
%       [  1 - plot when completed                                          ]
%       [  2 - plot incrementally at each iteration                         ]
%       [  3 - same as 2, with 'pause' at each iteration                    ]
%   cpf.plot.bus            <empty>     index of bus whose voltage is to be
%                                       plotted
%   cpf.user_callback       <empty>     string containing the name of a user
%                                       callback function, or struct with
%                                       function name, and optional priority
%                                       and/or args, or cell array of such
%                                       strings and/or structs, see
%                                       'help cpf_default_callback' for details
%
%Optimal Power Flow options:
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%   opf.ac.solver           'DEFAULT'   AC optimal power flow solver
%       [ 'DEFAULT' - choose solver based on availability in the following  ]
%       [             order: 'PDIPM', 'MIPS'                                ]
%       [ 'MIPS'    - MIPS, MATPOWER Interior Point Solver, primal/dual     ]
%       [             interior point method (pure Matlab)                   ]
%       [ 'FMINCON' - MATLAB Optimization Toolbox, FMINCON                  ]
%       [ 'IPOPT'   - IPOPT, requires MEX interface to IPOPT solver         ]
%       [             available from:                                       ]
%       [                 http://www.coin-or.org/projects/Ipopt.xml         ]
%       [ 'KNITRO'  - KNITRO, requires MATLAB Optimization Toolbox and      ]
%       [             KNITRO libraries available from: http://www.ziena.com/]
%       [ 'MINOPF'  - MINOPF, MINOS-based solver, requires optional         ]
%       [             MEX-based MINOPF package, available from:             ]
%       [                   http://www.pserc.cornell.edu/minopf/            ]
%       [ 'PDIPM'   - PDIPM, primal/dual interior point method, requires    ]
%       [             optional MEX-based TSPOPF package, available from:    ]
%       [                   http://www.pserc.cornell.edu/tspopf/            ]
%       [ 'SDPOPF'  - SDPOPF, solver based on semidefinite relaxation of    ]
%       [             OPF problem, requires optional packages:              ]
%       [               SDP_PF, available in extras/sdp_pf                  ]
%       [               YALMIP, available from:                             ]
%       [                   http://users.isy.liu.se/johanl/yalmip/          ]
%       [               SDP solver such as SeDuMi, available from:          ]
%       [                   http://sedumi.ie.lehigh.edu/                    ]
%       [ 'TRALM'   - TRALM, trust region based augmented Langrangian       ]
%       [             method, requires TSPOPF (see 'PDIPM')                 ]
%   opf.dc.solver           'DEFAULT'   DC optimal power flow solver
%       [ 'DEFAULT' - choose solver based on availability in the following  ]
%       [             order: 'GUROBI', 'CPLEX', 'MOSEK', 'OT',              ]
%       [             'GLPK' (linear costs only), 'BPMPD', 'MIPS'           ]
%       [ 'MIPS'    - MIPS, MATPOWER Interior Point Solver, primal/dual     ]
%       [             interior point method (pure Matlab)                   ]
%       [ 'BPMPD'   - BPMPD, requires optional MEX-based BPMPD_MEX package  ]
%       [             available from: http://www.pserc.cornell.edu/bpmpd/   ]
%       [ 'CLP'     - CLP, requires interface to COIN-OP LP solver          ]
%       [             available from:http://www.coin-or.org/projects/Clp.xml]
%       [ 'CPLEX'   - CPLEX, requires CPLEX solver available from:          ]
%       [             http://www.ibm.com/software/integration/ ...          ]
%       [                                 ... optimization/cplex-optimizer/ ]
%       [ 'GLPK'    - GLPK, requires interface to GLPK solver               ]
%       [             available from: http://www.gnu.org/software/glpk/     ]
%       [             (GLPK does not work with quadratic cost functions)    ]
%       [ 'GUROBI'  - GUROBI, requires Gurobi optimizer (v. 5+)             ]
%       [             available from: http://www.gurobi.com/                ]
%       [ 'IPOPT'   - IPOPT, requires MEX interface to IPOPT solver         ]
%       [             available from:                                       ]
%       [                 http://www.coin-or.org/projects/Ipopt.xml         ]
%       [ 'MOSEK'   - MOSEK, requires Matlab interface to MOSEK solver      ]
%       [             available from: http://www.mosek.com/                 ]
%       [ 'OT'      - MATLAB Optimization Toolbox, QUADPROG, LINPROG        ]
%   opf.violation           5e-6        constraint violation tolerance
%   opf.use_vg              0           respect gen voltage setpt     [ 0-1 ]
%       [ 0 - use specified bus Vmin & Vmax, and ignore gen Vg              ]
%       [ 1 - replace specified bus Vmin & Vmax by corresponding gen Vg     ]
%       [ between 0 and 1 - use a weighted average of the 2 options         ]
%   opf.flow_lim            'S'         quantity limited by branch flow
%                                       constraints
%       [ 'S' - apparent power flow (limit in MVA)                          ]
%       [ 'P' - active power flow (limit in MW)                             ]
%       [ 'I' - current magnitude (limit in MVA at 1 p.u. voltage)          ]
%   opf.ignore_angle_lim    0           angle diff limits for branches
%       [ 0 - include angle difference limits, if specified                 ]
%       [ 1 - ignore angle difference limits even if specified              ]
%   opf.init_from_mpc       -1          specify whether to use current state
%                                       in MATPOWER case to initialize OPF
%                                       (currently supported only for Ipopt,
%                                        Knitro and MIPS solvers)
%       [  -1 - MATPOWER decides, based on solver/algorithm                 ]
%       [   0 - ignore current state when initializing OPF                  ]
%       [   1 - use current state to initialize OPF                         ]
%   opf.return_raw_der      0           for AC OPF, return constraint and
%                                       derivative info in results.raw
%                                       (in fields g, dg, df, d2f) [ 0 or 1 ]
%
%Output options:
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%   verbose                 1           amount of progress info printed
%       [   0 - print no progress info                                      ]
%       [   1 - print a little progress info                                ]
%       [   2 - print a lot of progress info                                ]
%       [   3 - print all progress info                                     ]
%   out.all                 -1          controls pretty-printing of results
%       [  -1 - individual flags control what prints                        ]
%       [   0 - do not print anything (overrides individual flags, ignored  ]
%       [       for files specified as FNAME arg to runpf(), runopf(), etc.)]
%       [   1 - print everything (overrides individual flags)               ]
%   out.sys_sum             1           print system summary       [ 0 or 1 ]
%   out.area_sum            0           print area summaries       [ 0 or 1 ]
%   out.bus                 1           print bus detail           [ 0 or 1 ]
%   out.branch              1           print branch detail        [ 0 or 1 ]
%   out.gen                 0           print generator detail     [ 0 or 1 ]
%   out.lim.all             -1          controls constraint info output
%       [  -1 - individual flags control what constraint info prints        ]
%       [   0 - no constraint info (overrides individual flags)             ]
%       [   1 - binding constraint info (overrides individual flags)        ]
%       [   2 - all constraint info (overrides individual flags)            ]
%   out.lim.v               1           control voltage limit info
%       [   0 - do not print                                                ]
%       [   1 - print binding constraints only                              ]
%       [   2 - print all constraints                                       ]
%       [   (same options for OUT_LINE_LIM, OUT_PG_LIM, OUT_QG_LIM)         ]
%   out.lim.line            1           control line flow limit info
%   out.lim.pg              1           control gen active power limit info
%   out.lim.qg              1           control gen reactive pwr limit info
%   out.force               0           print results even if success
%                                       flag = 0                   [ 0 or 1 ]
%   out.suppress_detail     -1          suppress all output but system summary
%       [  -1 - suppress details for large systems (> 500 buses)            ]
%       [   0 - do not suppress any output specified by other flags         ]
%       [   1 - suppress all output except system summary section           ]
%       [       (overrides individual flags, but not out.all = 1)           ]
%
%Solver specific options:
%       name                    default     description [options]
%   -----------------------    ---------   ----------------------------------
%   MIPS:
%       mips.linsolver          ''          linear system solver
%           [   '' or '\'   build-in backslash \ operator (e.g. x = A \ b)  ]
%           [   'PARDISO'   PARDISO solver (if available)                   ]
%       mips.feastol            0           feasibility (equality) tolerance
%                                           (set to opf.violation by default)
%       mips.gradtol            1e-6        gradient tolerance
%       mips.comptol            1e-6        complementary condition
%                                           (inequality) tolerance
%       mips.costtol            1e-6        optimality tolerance
%       mips.max_it             150         maximum number of iterations
%       mips.step_control       0           enable step-size cntrl [ 0 or 1 ]
%       mips.sc.red_it          20          maximum number of reductions per
%                                           iteration with step control
%       mips.xi                 0.99995     constant used in alpha updates*
%       mips.sigma              0.1         centering parameter*
%       mips.z0                 1           used to initialize slack variables*
%       mips.alpha_min          1e-8        returns "Numerically Failed" if
%                                           either alpha parameter becomes
%                                           smaller than this value*
%       mips.rho_min            0.95        lower bound on rho_t*
%       mips.rho_max            1.05        upper bound on rho_t*
%       mips.mu_threshold       1e-5        KT multipliers smaller than this
%                                           value for non-binding constraints
%                                           are forced to zero
%       mips.max_stepsize       1e10        returns "Numerically Failed" if the
%                                           2-norm of the reduced Newton step
%                                           exceeds this value*
%           * See the corresponding Appendix in the manual for details.
%
%   CPLEX:
%       cplex.lpmethod          0           solution algorithm for LP problems
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex                                          ]
%           [   2 - dual simplex                                            ]
%           [   3 - network simplex                                         ]
%           [   4 - barrier                                                 ]
%           [   5 - sifting                                                 ]
%           [   6 - concurrent (dual, barrier, and primal)                  ]
%       cplex.qpmethod          0           solution algorithm for QP problems
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex optimizer                                ]
%           [   2 - dual simplex optimizer                                  ]
%           [   3 - network optimizer                                       ]
%           [   4 - barrier optimizer                                       ]
%       cplex.opts              <empty>     see CPLEX_OPTIONS for details
%       cplex.opt_fname         <empty>     see CPLEX_OPTIONS for details
%       cplex.opt               0           see CPLEX_OPTIONS for details
%
%   FMINCON:
%       fmincon.alg             4           algorithm used by fmincon() for OPF
%                                           for Opt Toolbox 4 and later
%            [  1 - active-set (not suitable for large problems)            ]
%            [  2 - interior-point, w/default 'bfgs' Hessian approx         ]
%            [  3 - interior-point, w/ 'lbfgs' Hessian approx               ]
%            [  4 - interior-point, w/exact user-supplied Hessian           ]
%            [  5 - interior-point, w/Hessian via finite differences        ]
%            [  6 - sqp (not suitable for large problems)                   ]
%       fmincon.tol_x           1e-4        termination tol on x
%       fmincon.tol_f           1e-4        termination tol on f
%       fmincon.max_it          0           maximum number of iterations
%                                                           [  0 => default ]
%
%   GUROBI:
%       gurobi.method           0           solution algorithm (Method)
%           [  -1 - automatic, let Gurobi decide                            ]
%           [   0 - primal simplex                                          ]
%           [   1 - dual simplex                                            ]
%           [   2 - barrier                                                 ]
%           [   3 - concurrent (LP only)                                    ]
%           [   4 - deterministic concurrent (LP only)                      ]
%       gurobi.timelimit        Inf         maximum time allowed (TimeLimit)
%       gurobi.threads          0           max number of threads (Threads)
%       gurobi.opts             <empty>     see GUROBI_OPTIONS for details
%       gurobi.opt_fname        <empty>     see GUROBI_OPTIONS for details
%       gurobi.opt              0           see GUROBI_OPTIONS for details
%
%   IPOPT:
%       ipopt.opts              <empty>     see IPOPT_OPTIONS for details
%       ipopt.opt_fname         <empty>     see IPOPT_OPTIONS for details
%       ipopt.opt               0           see IPOPT_OPTIONS for details
%
%   KNITRO:
%       knitro.tol_x            1e-4        termination tol on x
%       knitro.tol_f            1e-4        termination tol on f
%       knitro.opt_fname        <empty>     name of user-supplied native
%                                           KNITRO options file that overrides
%                                           all other options
%       knitro.opt              0           if knitro.opt_fname is empty and
%                                           knitro.opt is a non-zero integer N
%                                           then knitro.opt_fname is auto-
%                                           generated as:
%                                           'knitro_user_options_N.txt'
%
%   LINPROG:
%       linprog                 <empty>     LINPROG options passed to
%                                           OPTIMOPTIONS or OPTIMSET.
%                                           see LINPROG in the Optimization
%                                           Toolbox for details
%
%   MINOPF:
%       minopf.feastol          0 (1e-3)    primal feasibility tolerance
%                                           (set to opf.violation by default)
%       minopf.rowtol           0 (1e-3)    row tolerance
%       minopf.xtol             0 (1e-4)    x tolerance
%       minopf.majdamp          0 (0.5)     major damping parameter
%       minopf.mindamp          0 (2.0)     minor damping parameter
%       minopf.penalty          0 (1.0)     penalty parameter
%       minopf.major_it         0 (200)     major iterations
%       minopf.minor_it         0 (2500)    minor iterations
%       minopf.max_it           0 (2500)    iterations limit
%       minopf.verbosity        -1          amount of progress info printed
%           [  -1 - controlled by 'verbose' option                          ]
%           [   0 - print nothing                                           ]
%           [   1 - print only termination status message                   ]
%           [   2 - print termination status and screen progress            ]
%           [   3 - print screen progress, report file (usually fort.9)     ]
%       minopf.core             0 (1200*nb + 2*(nb+ng)^2) memory allocation
%       minopf.supbasic_lim     0 (2*nb + 2*ng) superbasics limit
%       minopf.mult_price       0 (30)      multiple price
%
%   MOSEK:
%       mosek.lp_alg            0           solution algorithm
%                                               (MSK_IPAR_OPTIMIZER)
%           for MOSEK 8.x ...        (see MOSEK_SYMBCON for a "better way")
%           [   0 - automatic: let MOSEK choose                             ]
%           [   1 - dual simplex                                            ]
%           [   2 - automatic: let MOSEK choose                             ]
%           [   3 - automatic simplex (MOSEK chooses which simplex method)  ]
%           [   4 - interior point                                          ]
%           [   6 - primal simplex                                          ]
%       mosek.max_it            0 (400)     interior point max iterations
%                                               (MSK_IPAR_INTPNT_MAX_ITERATIONS)
%       mosek.gap_tol           0 (1e-8)    interior point relative gap tol
%                                               (MSK_DPAR_INTPNT_TOL_REL_GAP)
%       mosek.max_time          0 (-1)      maximum time allowed
%                                               (MSK_DPAR_OPTIMIZER_MAX_TIME)
%       mosek.num_threads       0 (1)       max number of threads
%                                               (MSK_IPAR_INTPNT_NUM_THREADS)
%       mosek.opts              <empty>     see MOSEK_OPTIONS for details
%       mosek.opt_fname         <empty>     see MOSEK_OPTIONS for details
%       mosek.opt               0           see MOSEK_OPTIONS for details
%
%   QUADPROG:
%       quadprog                <empty>     QUADPROG options passed to
%                                           OPTIMOPTIONS or OPTIMSET.
%                                           see QUADPROG in the Optimization
%                                           Toolbox for details
%
%   TSPOPF:
%       pdipm.feastol           0           feasibility (equality) tolerance
%                                           (set to opf.violation by default)
%       pdipm.gradtol           1e-6        gradient tolerance
%       pdipm.comptol           1e-6        complementary condition
%                                           (inequality) tolerance
%       pdipm.costtol           1e-6        optimality tolerance
%       pdipm.max_it            150         maximum number of iterations
%       pdipm.step_control      0           enable step-size cntrl [ 0 or 1 ]
%       pdipm.sc.red_it         20          maximum number of reductions per
%                                           iteration with step control
%       pdipm.sc.smooth_ratio   0.04        piecewise linear curve smoothing
%                                           ratio
%
%       tralm.feastol           0           feasibility tolerance
%                                           (set to opf.violation by default)
%       tralm.primaltol         5e-4        primal variable tolerance
%       tralm.dualtol           5e-4        dual variable tolerance
%       tralm.costtol           1e-5        optimality tolerance
%       tralm.major_it          40          maximum number of major iterations
%       tralm.minor_it          40          maximum number of minor iterations
%       tralm.smooth_ratio      0.04        piecewise linear curve smoothing
%                                           ratio
%
%Experimental Options:
%   exp.sys_wide_zip_loads.pw   <empty>     1 x 3 vector of active load fraction
%                                           to be modeled as constant power,
%                                           constant current and constant
%                                           impedance, respectively, where
%                                           <empty> means use [1 0 0]
%   exp.sys_wide_zip_loads.qw   <empty>     same for reactive power, where
%                                           <empty> means use same value as
%                                           for 'pw'

%   MATPOWER
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% some constants
N = 124;                %% number of options in old-style vector (MATPOWER 4.1)
N40 = 116;              %% dimension of MATPOWER 4.0 options vector
N32 = 93;               %% dimension of MATPOWER 3.2 options vector
v = mpoption_version;   %% version number of MATPOWER options struct

%% initialize flags and arg counter
have_opt0 = 0;          %% existing options struct or vector provided?
have_old_style_ov = 0;  %% override options using old-style names?
return_old_style = 0;   %% return value as old-style vector?
k = 1;
if nargin > 0
    opt0 = varargin{k};
    if isstruct(opt0) && isfield(opt0, 'v')         %% options struct
        have_opt0 = 1;
        k = k + 1;
    elseif isnumeric(opt0) && size(opt0, 2) == 1    %% options vector
        nn = size(opt0, 1);
        if ismember(nn, [N N40 N32])                %% of valid size
            %% expand smaller option vectors (from before MATPOWER 4.1)
            if nn < N
                optv = mpoption_old();
                opt0(nn+1:N) = optv(nn+1:N);
            end
            have_opt0 = 1;
            k = k + 1;
        end
    end
end

%% create base options vector to which overrides are made
if have_opt0
    if isstruct(opt0)               %% it's already a valid options struct
        if DEBUG, fprintf('OPT0 is a valid options struct\n'); end
        if opt0.v < v
            %% convert older version to current version
            opt_d = mpoption_default();
            if opt0.v == 1          %% convert version 1 to 2
                if isfield(opt_d, 'linprog')
                    opt0.lingprog = opt_d.linprog;
                end
                if isfield(opt_d, 'quadprog')
                    opt0.quadprog = opt_d.quadprog;
                end
            end
            if opt0.v <= 2          %% convert version 2 to 3
                opt0.out.suppress_detail = opt_d.out.suppress_detail;
            end
            %if opt0.v <= 3          %% convert version 3 to 4
                %% new mips options were all optional, no conversion needed
            %end
            if opt0.v <= 4          %% convert version 4 to 5
                opt0.opf.init_from_mpc = opt_d.opf.init_from_mpc;
            end
            if opt0.v <= 5          %% convert version 5 to 6
                if isfield(opt_d, 'clp')
                    opt0.clp = opt_d.clp;
                end
            end
            if opt0.v <= 6          %% convert version 6 to 7
                if isfield(opt_d, 'intlinprog')
                    opt0.intlinprog = opt_d.intlinprog;
                end
            end
            if opt0.v <= 7          %% convert version 7 to 8
                opt0.mips.linsolver = opt_d.mips.linsolver;
            end
            if opt0.v <= 8          %% convert version 8 to 9
                opt0.exp.sys_wide_zip_loads = opt_d.exp.sys_wide_zip_loads;
            end
            if opt0.v <= 9          %% convert version 9 to 10
                opt0.most = opt_d.most;
            end
            if opt0.v <= 10         %% convert version 10 to 11
                opt0.cpf.enforce_p_lims     = opt_d.cpf.enforce_p_lims;
                opt0.cpf.enforce_q_lims     = opt_d.cpf.enforce_q_lims;
                opt0.cpf.adapt_step_damping = opt_d.cpf.adapt_step_damping;
                opt0.cpf.target_lam_tol     = opt_d.cpf.target_lam_tol;
                opt0.cpf.nose_tol           = opt_d.cpf.nose_tol;
                opt0.cpf.p_lims_tol         = opt_d.cpf.p_lims_tol;
                opt0.cpf.q_lims_tol         = opt_d.cpf.q_lims_tol;
                if (~isempty(opt0.cpf.user_callback_args) && ...
                        ~isstruct(opt0.cpf.user_callback_args)) || ...
                        (isstruct(opt0.cpf.user_callback_args) && ...
                         ~isempty(fields(opt0.cpf.user_callback_args)))
                    warning('The ''cpf.user_callback_args'' option has been removed. Please include the args in a struct in ''cpf.user_callback'' instead.')
                end
                opt0.cpf = rmfield(opt0.cpf, 'user_callback_args');
            end
            if opt0.v <= 11         %% convert version 11 to 12
                opt0.cpf.use_vg = opt_d.cpf.use_vg;
            end
            opt0.v = v;
        end
        opt = opt0;
    else                            %% convert from old-style options vector
        if DEBUG, fprintf('OPT0 is a old-style options vector\n'); end
        opt = mpoption_v2s(opt0);
    end
else                                %% use default options struct as base
    if DEBUG, fprintf('no OPT0, starting with default options struct\n'); end
    opt = mpoption_default();
end


%% do we have OVERRIDES or NAME/VALUE pairs
ov = [];
if nargin - k == 0          %% looking at last arg, must be OVERRIDES
    if isstruct(varargin{k})        %% OVERRIDES provided as struct
        if DEBUG, fprintf('OVERRIDES struct\n'); end
        ov = varargin{k};
    elseif ischar(varargin{k})      %% OVERRIDES provided as file/function name
        if DEBUG, fprintf('OVERRIDES file/function name\n'); end
        try
            ov = feval(varargin{k});
        catch
            error('mpoption: Unable to load MATPOWER options from ''%s''', varargin{k});
        end
        if ~isstruct(ov)
            error('mpoption: calling ''%s'' did not return a struct', varargin{k});
        end
    elseif isempty(varargin{k})
        return_old_style = 1;
    else
        error('mpoption: OVERRIDES must be a struct or the name of a function that returns a struct');
    end
elseif nargin - k > 0 && mod(nargin-k, 2)   %% even number of remaining args
    if DEBUG, fprintf('NAME/VALUE pairs override defaults\n'); end
    %% process NAME/VALUE pairs
    if strcmp(varargin{k}, upper(varargin{k}))  %% old-style, all UPPERCASE option pairs
        %% NOTE: new top-level option fields cannot be all uppercase
        if ~have_opt0
            opt_v = mpoption_old(varargin{:});  %% create modified vector ...
            opt = mpoption_v2s(opt_v);          %% ... then convert
        else
            have_old_style_ov = 1;
            %% convert pairs to struct
            while k < nargin
                name = varargin{k};
                val  = varargin{k+1};
                k = k + 2;
                ov.(name) = val;
            end
        end
    else                                        %% new option pairs
        %% convert pairs to struct
        while k < nargin
            name = varargin{k};
            val  = varargin{k+1};
            k = k + 2;
            c = regexp(name, '([^\.]*)', 'tokens');
            s = struct();
            for i = 1:length(c)
                s(i).type = '.';
                s(i).subs = c{i}{1};
            end
            ov = subsasgn(ov, s, val);
        end
    end
elseif nargin == 0 || nargin == 1
    if DEBUG, fprintf('no OVERRIDES, return default options struct or converted OPT0 vector\n'); end
else
    error('mpoption: invalid calling syntax, see ''help mpoption'' to double-check the valid options');
end

%% apply overrides
if ~isempty(ov)
    if have_old_style_ov
        opt = apply_old_mpoption_overrides(opt, ov);
    else
        persistent nsc_opt;     %% cache this to speed things up
        if ~isstruct(nsc_opt)
            vf = nested_struct_copy(mpoption_default(), mpoption_info_mips('V'));
            vf = nested_struct_copy(vf, mpoption_optional_fields());
            ex = struct(...
                'name', {}, ...
                'check', {}, ...
                'copy_mode', {} ...
            );
            %% add exceptions for optional packages
            opt_pkgs = mpoption_optional_pkgs();
            n = length(ex);
            for k = 1:length(opt_pkgs)
                fname = ['mpoption_info_' opt_pkgs{k}];
                if exist(fname, 'file') == 2
                    opt_ex = feval(fname, 'E');
                    nex = length(opt_ex);
                    if ~isempty(opt_ex)
                        for j = 1:nex
                            ex(n+j).name = opt_ex(j).name;
                        end
                        if isfield(opt_ex, 'check')
                            for j = 1:nex
                                ex(n+j).check = opt_ex(j).check;
                            end
                        end
                        if isfield(opt_ex, 'copy_mode')
                            for j = 1:nex
                                ex(n+j).copy_mode = opt_ex(j).copy_mode;
                            end
                        end
                        if isfield(opt_ex, 'valid_fields')
                            for j = 1:nex
                                ex(n+j).valid_fields = opt_ex(j).valid_fields;
                            end
                        end
                        n = n + nex;
                    end
                end
            end
            nsc_opt = struct('check', 1, 'valid_fields', vf, 'exceptions', ex);
        end
%         if have_fcn('catchme')
%             try
%                 opt = nested_struct_copy(opt, ov, nsc_opt);
%             catch me
%                 str = strrep(me.message, 'field', 'option');
%                 str = strrep(str, 'nested_struct_copy', 'mpoption');
%                 error(str);
%             end
%         else
            try
                opt = nested_struct_copy(opt, ov, nsc_opt);
            catch
                me = lasterr;
                str = strrep(me, 'field', 'option');
                str = strrep(str, 'nested_struct_copy', 'mpoption');
                error(str);
            end
%         end
    end
end
if return_old_style
    opt = mpoption_s2v(opt);
end


%%-------------------------------------------------------------------
function opt = apply_old_mpoption_overrides(opt0, ov)
%
%   OPT0 is assumed to already have all of the fields and sub-fields found
%   in the default options struct.

%% initialize output
opt = opt0;

errstr = 'mpoption: %g is not a valid value for the old-style ''%s'' option';
fields = fieldnames(ov);
for f = 1:length(fields)
    ff = fields{f};
    switch ff
        case 'PF_ALG'
            switch ov.(ff)
                case 1
                    opt.pf.alg = 'NR';      %% Newton's method
                case 2
                    opt.pf.alg = 'FDXB';    %% fast-decoupled (XB version)
                case 3
                    opt.pf.alg = 'FDBX';    %% fast-decoupled (BX version)
                case 4
                    opt.pf.alg = 'GS';      %% Gauss-Seidel
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'PF_TOL'
            opt.pf.tol = ov.(ff);
        case 'PF_MAX_IT'
            opt.pf.nr.max_it = ov.(ff);
        case 'PF_MAX_IT_FD'
            opt.pf.fd.max_it = ov.(ff);
        case 'PF_MAX_IT_GS'
            opt.pf.gs.max_it = ov.(ff);
        case 'ENFORCE_Q_LIMS'
            opt.pf.enforce_q_lims = ov.(ff);
        case 'PF_DC'
            switch ov.(ff)
                case 0
                    opt.model = 'AC';
                case 1
                    opt.model = 'DC';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_ALG'
            switch ov.(ff)
                case 0
                    opt.opf.ac.solver = 'DEFAULT';
                case 500
                    opt.opf.ac.solver = 'MINOPF';
                case 520
                    opt.opf.ac.solver = 'FMINCON';
                case {540, 545}
                    opt.opf.ac.solver = 'PDIPM';
                    if ov.(ff) == 545
                        opt.pdipm.step_control = 1;
                    else
                        opt.pdipm.step_control = 0;
                    end
                case 550
                    opt.opf.ac.solver = 'TRALM';
                case {560, 565}
                    opt.opf.ac.solver = 'MIPS';
                    if ov.(ff) == 565
                        opt.mips.step_control = 1;
                    else
                        opt.mips.step_control = 0;
                    end
                case 580
                    opt.opf.ac.solver = 'IPOPT';
                case 600
                    opt.opf.ac.solver = 'KNITRO';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_VIOLATION'
            opt.opf.violation = ov.(ff);
        case 'CONSTR_TOL_X'
            opt.fmincon.tol_x = ov.(ff);
            opt.knitro.tol_x = ov.(ff);
        case 'CONSTR_TOL_F'
            opt.fmincon.tol_f = ov.(ff);
            opt.knitro.tol_f = ov.(ff);
        case 'CONSTR_MAX_IT'
            opt.fmincon.max_it = ov.(ff);
        case 'OPF_FLOW_LIM'
            switch ov.(ff)
                case 0
                    opt.opf.flow_lim = 'S';   %% apparent power (MVA)
                case 1
                    opt.opf.flow_lim = 'P';   %% real power (MW)
                case 2
                    opt.opf.flow_lim = 'I';   %% current magnitude (MVA @ 1 p.u. voltage)
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_IGNORE_ANG_LIM'
            opt.opf.ignore_angle_lim = ov.(ff);
        case 'OPF_ALG_DC'
            switch ov.(ff)
                case 0
                    opt.opf.dc.solver = 'DEFAULT';
                case 100
                    opt.opf.dc.solver = 'BPMPD';
                case {200, 250}
                    opt.opf.dc.solver = 'MIPS';
                    if ov.(ff) == 250
                        opt.mips.step_control = 1;
                    else
                        opt.mips.step_control = 0;
                    end
                case 300
                    opt.opf.dc.solver = 'OT';     %% QUADPROG, LINPROG
                case 400
                    opt.opf.dc.solver = 'IPOPT';
                case 500
                    opt.opf.dc.solver = 'CPLEX';
                case 600
                    opt.opf.dc.solver = 'MOSEK';
                case 700
                    opt.opf.dc.solver = 'GUROBI';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'VERBOSE'
            opt.verbose = ov.(ff);
        case 'OUT_ALL'
            opt.out.all = ov.(ff);
        case 'OUT_SYS_SUM'
            opt.out.sys_sum = ov.(ff);
        case 'OUT_AREA_SUM'
            opt.out.area_sum = ov.(ff);
        case 'OUT_BUS'
            opt.out.bus = ov.(ff);
        case 'OUT_BRANCH'
            opt.out.branch = ov.(ff);
        case 'OUT_GEN'
            opt.out.gen = ov.(ff);
        case 'OUT_ALL_LIM'
            opt.out.lim.all = ov.(ff);
        case 'OUT_V_LIM'
            opt.out.lim.v = ov.(ff);
        case 'OUT_LINE_LIM'
            opt.out.lim.line = ov.(ff);
        case 'OUT_PG_LIM'
            opt.out.lim.pg = ov.(ff);
        case 'OUT_QG_LIM'
            opt.out.lim.qg = ov.(ff);
        case 'OUT_FORCE'
            opt.out.force = ov.(ff);
        case 'RETURN_RAW_DER'
            opt.opf.return_raw_der = ov.(ff);
        case 'FMC_ALG'
            opt.fmincon.alg = ov.(ff);
        case 'KNITRO_OPT'
            opt.knitro.opt = ov.(ff);
        case 'IPOPT_OPT'
            opt.ipopt.opt = ov.(ff);
        case 'MNS_FEASTOL'
            opt.minopf.feastol = ov.(ff);
        case 'MNS_ROWTOL'
            opt.minopf.rowtol = ov.(ff);
        case 'MNS_XTOL'
            opt.minopf.xtol = ov.(ff);
        case 'MNS_MAJDAMP'
            opt.minopf.majdamp = ov.(ff);
        case 'MNS_MINDAMP'
            opt.minopf.mindamp = ov.(ff);
        case 'MNS_PENALTY_PARM'
            opt.minopf.penalty = ov.(ff);
        case 'MNS_MAJOR_IT'
            opt.minopf.major_it = ov.(ff);
        case 'MNS_MINOR_IT'
            opt.minopf.minor_it = ov.(ff);
        case 'MNS_MAX_IT'
            opt.minopf.max_it = ov.(ff);
        case 'MNS_VERBOSITY'
            opt.minopf.verbosity = ov.(ff);
        case 'MNS_CORE'
            opt.minopf.core = ov.(ff);
        case 'MNS_SUPBASIC_LIM'
            opt.minopf.supbasic_lim = ov.(ff);
        case 'MNS_MULT_PRICE'
            opt.minopf.mult_price = ov.(ff);
        case 'FORCE_PC_EQ_P0'
            opt.sopf.force_Pc_eq_P0 = ov.(ff);
        case 'PDIPM_FEASTOL'
            opt.mips.feastol = ov.(ff);
            opt.pdipm.feastol = ov.(ff);
        case 'PDIPM_GRADTOL'
            opt.mips.gradtol = ov.(ff);
            opt.pdipm.gradtol = ov.(ff);
        case 'PDIPM_COMPTOL'
            opt.mips.comptol = ov.(ff);
            opt.pdipm.comptol = ov.(ff);
        case 'PDIPM_COSTTOL'
            opt.mips.costtol = ov.(ff);
            opt.pdipm.costtol = ov.(ff);
        case 'PDIPM_MAX_IT'
            opt.mips.max_it = ov.(ff);
            opt.pdipm.max_it = ov.(ff);
        case 'SCPDIPM_RED_IT'
            opt.mips.sc.red_it = ov.(ff);
            opt.pdipm.sc.red_it = ov.(ff);
        case 'TRALM_FEASTOL'
            opt.tralm.feastol = ov.(ff);
        case 'TRALM_PRIMETOL'
            opt.tralm.primaltol = ov.(ff);
        case 'TRALM_DUALTOL'
            opt.tralm.dualtol = ov.(ff);
        case 'TRALM_COSTTOL'
            opt.tralm.costtol = ov.(ff);
        case 'TRALM_MAJOR_IT'
            opt.tralm.major_it = ov.(ff);
        case 'TRALM_MINOR_IT'
            opt.tralm.minor_it = ov.(ff);
        case 'SMOOTHING_RATIO'
            opt.pdipm.sc.smooth_ratio = ov.(ff);
            opt.tralm.smooth_ratio = ov.(ff);
        case 'CPLEX_LPMETHOD'
            opt.cplex.lpmethod = ov.(ff);
        case 'CPLEX_QPMETHOD'
            opt.cplex.qpmethod = ov.(ff);
        case 'CPLEX_OPT'
            opt.cplex.opt = ov.(ff);
        case 'MOSEK_LP_ALG'
            opt.mosek.lp_alg = ov.(ff);
        case 'MOSEK_MAX_IT'
            opt.mosek.max_it = ov.(ff);
        case 'MOSEK_GAP_TOL'
            opt.mosek.gap_tol = ov.(ff);
        case 'MOSEK_MAX_TIME'
            opt.mosek.max_time = ov.(ff);
        case 'MOSEK_NUM_THREADS'
            opt.mosek.num_threads = ov.(ff);
        case 'MOSEK_OPT'
            opt.mosek.opt = ov.(ff);
        case 'GRB_METHOD'
            opt.gurobi.method = ov.(ff);
        case 'GRB_TIMELIMIT'
            opt.gurobi.timelimit = ov.(ff);
        case 'GRB_THREADS'
            opt.gurobi.threads = ov.(ff);
        case 'GRB_OPT'
            opt.gurobi.opt = ov.(ff);
        otherwise
            error('mpoption: ''%s'' is not a valid old-style option name', ff);
    end
end
% ov


%%-------------------------------------------------------------------
function opt_s = mpoption_v2s(opt_v)
if DEBUG, fprintf('mpoption_v2s()\n'); end
opt_s = mpoption_default();
errstr = 'mpoption: %g is not a valid value for the old-style ''%s'' option';
switch opt_v(1)                                 %% PF_ALG
    case 1
        opt_s.pf.alg = 'NR';        %% Newton's method
    case 2
        opt_s.pf.alg = 'FDXB';      %% fast-decoupled (XB version)
    case 3
        opt_s.pf.alg = 'FDBX';      %% fast-decoupled (BX version)
    case 4
        opt_s.pf.alg = 'GS';        %% Gauss-Seidel
    otherwise
        error(errstr, opt_v(1), 'PF_ALG');
end
opt_s.pf.tol                = opt_v(2);         %% PF_TOL
opt_s.pf.nr.max_it          = opt_v(3);         %% PF_MAX_IT
opt_s.pf.fd.max_it          = opt_v(4);         %% PF_MAX_IT_FD
opt_s.pf.gs.max_it          = opt_v(5);         %% PF_MAX_IT_GS
opt_s.pf.enforce_q_lims     = opt_v(6);         %% ENFORCE_Q_LIMS
switch opt_v(10)                                %% PF_DC
    case 0
        opt_s.model = 'AC';
    case 1
        opt_s.model = 'DC';
    otherwise
        error(errstr, opt_v(10), 'PF_DC');
end
switch opt_v(11)                                %% OPF_ALG
    case 0
        opt_s.opf.ac.solver = 'DEFAULT';
    case 500
        opt_s.opf.ac.solver = 'MINOPF';
    case 520
        opt_s.opf.ac.solver = 'FMINCON';
    case {540, 545}
        opt_s.opf.ac.solver = 'PDIPM';
    case 550
        opt_s.opf.ac.solver = 'TRALM';
    case {560, 565}
        opt_s.opf.ac.solver = 'MIPS';
    case 580
        opt_s.opf.ac.solver = 'IPOPT';
    case 600
        opt_s.opf.ac.solver = 'KNITRO';
    otherwise
        error(errstr, opt_v(11), 'OPF_ALG');
end
opt_s.opf.violation         = opt_v(16);        %% OPF_VIOLATION

opt_s.fmincon.tol_x         = opt_v(17);        %% CONSTR_TOL_X
opt_s.fmincon.tol_f         = opt_v(18);        %% CONSTR_TOL_F
opt_s.fmincon.max_it        = opt_v(19);        %% CONSTR_MAX_IT

opt_s.knitro.tol_x          = opt_v(17);        %% CONSTR_TOL_X
opt_s.knitro.tol_f          = opt_v(18);        %% CONSTR_TOL_F

switch opt_v(24)                                %% OPF_FLOW_LIM
    case 0
        opt_s.opf.flow_lim = 'S';   %% apparent power (MVA)
    case 1
        opt_s.opf.flow_lim = 'P';   %% real power (MW)
    case 2
        opt_s.opf.flow_lim = 'I';   %% current magnitude (MVA @ 1 p.u. voltage)
    otherwise
        error(errstr, opt_v(10), 'PF_DC');
end

opt_s.opf.ignore_angle_lim  = opt_v(25);        %% OPF_IGNORE_ANG_LIM

switch opt_v(26)                                %% OPF_ALG_DC
    case 0
        opt_s.opf.dc.solver = 'DEFAULT';
    case 100
        opt_s.opf.dc.solver = 'BPMPD';
    case {200, 250}
        opt_s.opf.dc.solver = 'MIPS';
    case 300
        opt_s.opf.dc.solver = 'OT';     %% QUADPROG, LINPROG
    case 400
        opt_s.opf.dc.solver = 'IPOPT';
    case 500
        opt_s.opf.dc.solver = 'CPLEX';
    case 600
        opt_s.opf.dc.solver = 'MOSEK';
    case 700
        opt_s.opf.dc.solver = 'GUROBI';
    otherwise
        error(errstr, opt_v(26), 'OPF_ALG_DC');
end

opt_s.verbose               = opt_v(31);        %% VERBOSE
opt_s.out.all               = opt_v(32);        %% OUT_ALL
opt_s.out.sys_sum           = opt_v(33);        %% OUT_SYS_SUM
opt_s.out.area_sum          = opt_v(34);        %% OUT_AREA_SUM
opt_s.out.bus               = opt_v(35);        %% OUT_BUS
opt_s.out.branch            = opt_v(36);        %% OUT_BRANCH
opt_s.out.gen               = opt_v(37);        %% OUT_GEN
opt_s.out.lim.all           = opt_v(38);        %% OUT_ALL_LIM
opt_s.out.lim.v             = opt_v(39);        %% OUT_V_LIM
opt_s.out.lim.line          = opt_v(40);        %% OUT_LINE_LIM
opt_s.out.lim.pg            = opt_v(41);        %% OUT_PG_LIM
opt_s.out.lim.qg            = opt_v(42);        %% OUT_QG_LIM
opt_s.out.force             = opt_v(44);        %% OUT_FORCE

opt_s.opf.return_raw_der    = opt_v(52);        %% RETURN_RAW_DER

opt_s.fmincon.alg           = opt_v(55);        %% FMC_ALG
opt_s.knitro.opt            = opt_v(58);        %% KNITRO_OPT
opt_s.ipopt.opt             = opt_v(60);        %% IPOPT_OPT

opt_s.minopf.feastol        = opt_v(61);        %% MNS_FEASTOL
opt_s.minopf.rowtol         = opt_v(62);        %% MNS_ROWTOL
opt_s.minopf.xtol           = opt_v(63);        %% MNS_XTOL
opt_s.minopf.majdamp        = opt_v(64);        %% MNS_MAJDAMP
opt_s.minopf.mindamp        = opt_v(65);        %% MNS_MINDAMP
opt_s.minopf.penalty        = opt_v(66);        %% MNS_PENALTY_PARM
opt_s.minopf.major_it       = opt_v(67);        %% MNS_MAJOR_IT
opt_s.minopf.minor_it       = opt_v(68);        %% MNS_MINOR_IT
opt_s.minopf.max_it         = opt_v(69);        %% MNS_MAX_IT
opt_s.minopf.verbosity      = opt_v(70);        %% MNS_VERBOSITY
opt_s.minopf.core           = opt_v(71);        %% MNS_CORE
opt_s.minopf.supbasic_lim   = opt_v(72);        %% MNS_SUPBASIC_LIM
opt_s.minopf.mult_price     = opt_v(73);        %% MNS_MULT_PRICE

opt_s.sopf.force_Pc_eq_P0   = opt_v(80);        %% FORCE_PC_EQ_P0, for c3sopf

if (opt_v(11) == 565 && opt_v(10) == 0) || (opt_v(26) == 250 && opt_v(10) == 1)
    opt_s.mips.step_control = 1;
end
opt_s.mips.feastol          = opt_v(81);        %% PDIPM_FEASTOL
opt_s.mips.gradtol          = opt_v(82);        %% PDIPM_GRADTOL
opt_s.mips.comptol          = opt_v(83);        %% PDIPM_COMPTOL
opt_s.mips.costtol          = opt_v(84);        %% PDIPM_COSTTOL
opt_s.mips.max_it           = opt_v(85);        %% PDIPM_MAX_IT
opt_s.mips.sc.red_it        = opt_v(86);        %% SCPDIPM_RED_IT

opt_s.pdipm.feastol         = opt_v(81);        %% PDIPM_FEASTOL
opt_s.pdipm.gradtol         = opt_v(82);        %% PDIPM_GRADTOL
opt_s.pdipm.comptol         = opt_v(83);        %% PDIPM_COMPTOL
opt_s.pdipm.costtol         = opt_v(84);        %% PDIPM_COSTTOL
opt_s.pdipm.max_it          = opt_v(85);        %% PDIPM_MAX_IT
opt_s.pdipm.sc.red_it       = opt_v(86);        %% SCPDIPM_RED_IT
opt_s.pdipm.sc.smooth_ratio = opt_v(93);        %% SMOOTHING_RATIO
if opt_v(11) == 545 && opt_v(10) == 0
    opt_s.pdipm.step_control = 1;
end

opt_s.tralm.feastol         = opt_v(87);        %% TRALM_FEASTOL
opt_s.tralm.primaltol       = opt_v(88);        %% TRALM_PRIMETOL
opt_s.tralm.dualtol         = opt_v(89);        %% TRALM_DUALTOL
opt_s.tralm.costtol         = opt_v(90);        %% TRALM_COSTTOL
opt_s.tralm.major_it        = opt_v(91);        %% TRALM_MAJOR_IT
opt_s.tralm.minor_it        = opt_v(92);        %% TRALM_MINOR_IT
opt_s.tralm.smooth_ratio    = opt_v(93);        %% SMOOTHING_RATIO

opt_s.cplex.lpmethod        = opt_v(95);        %% CPLEX_LPMETHOD
opt_s.cplex.qpmethod        = opt_v(96);        %% CPLEX_QPMETHOD
opt_s.cplex.opt             = opt_v(97);        %% CPLEX_OPT

opt_s.mosek.lp_alg          = opt_v(111);       %% MOSEK_LP_ALG
opt_s.mosek.max_it          = opt_v(112);       %% MOSEK_MAX_IT
opt_s.mosek.gap_tol         = opt_v(113);       %% MOSEK_GAP_TOL
opt_s.mosek.max_time        = opt_v(114);       %% MOSEK_MAX_TIME
opt_s.mosek.num_threads     = opt_v(115);       %% MOSEK_NUM_THREADS
opt_s.mosek.opt             = opt_v(116);       %% MOSEK_OPT

opt_s.gurobi.method         = opt_v(121);       %% GRB_METHOD
opt_s.gurobi.timelimit      = opt_v(122);       %% GRB_TIMELIMIT
opt_s.gurobi.threads        = opt_v(123);       %% GRB_THREADS
opt_s.gurobi.opt            = opt_v(124);       %% GRB_OPT


%%-------------------------------------------------------------------
function opt_v = mpoption_s2v(opt_s)
if DEBUG, fprintf('mpoption_s2v()\n'); end
%% PF_ALG
old = mpoption_old;
switch upper(opt_s.pf.alg)
    case 'NR'
        PF_ALG = 1;
    case 'FDXB'
        PF_ALG = 2;
    case 'FDBX'
        PF_ALG = 3;
    case 'GS'
        PF_ALG = 4;
end

%% PF_DC
if strcmp(upper(opt_s.model), 'DC')
    PF_DC = 1;
else
    PF_DC = 0;
end

%% OPF_ALG
switch upper(opt_s.opf.ac.solver)
    case 'DEFAULT'
        OPF_ALG = 0;
    case 'MINOPF'
        OPF_ALG = 500;
    case 'FMINCON'
        OPF_ALG = 520;
    case 'PDIPM'
        if isfield(opt_s, 'pdipm') && opt_s.pdipm.step_control
            OPF_ALG = 545;
        else
            OPF_ALG = 540;
        end
    case 'TRALM'
        OPF_ALG = 550;
    case 'MIPS'
        if opt_s.mips.step_control
            OPF_ALG = 565;
        else
            OPF_ALG = 560;
        end
    case 'IPOPT'
        OPF_ALG = 580;
    case 'KNITRO'
        OPF_ALG = 600;
end

%% FMINCON, Knitro tol_x, tol_f, max_it
if strcmp(upper(opt_s.opf.ac.solver), 'KNITRO') && isfield(opt_s, 'knitro')
    CONSTR_TOL_X = opt_s.knitro.tol_x;
    CONSTR_TOL_F = opt_s.knitro.tol_f;
elseif isfield(opt_s, 'fmincon')
    CONSTR_TOL_X  = opt_s.fmincon.tol_x;
    CONSTR_TOL_F  = opt_s.fmincon.tol_f;
else
    CONSTR_TOL_X = old(17);
    CONSTR_TOL_F = old(18);
end
if isfield(opt_s, 'fmincon')
    CONSTR_MAX_IT   = opt_s.fmincon.max_it;
    FMC_ALG         = opt_s.fmincon.alg;
else
    CONSTR_MAX_IT   = old(19);
    FMC_ALG         = old(55);
end

%% OPF_FLOW_LIM
switch upper(opt_s.opf.flow_lim)
    case 'S'
        OPF_FLOW_LIM = 0;
    case 'P'
        OPF_FLOW_LIM = 1;
    case 'I'
        OPF_FLOW_LIM = 2;
end

%% OPF_ALG_DC
switch upper(opt_s.opf.dc.solver)
    case 'DEFAULT'
        OPF_ALG_DC = 0;
    case 'BPMPD'
        OPF_ALG_DC = 100;
    case 'MIPS'
        if opt_s.mips.step_control
            OPF_ALG_DC = 250;
        else
            OPF_ALG_DC = 200;
        end
    case 'OT'
        OPF_ALG_DC = 300;
    case 'IPOPT'
        OPF_ALG_DC = 400;
    case 'CPLEX'
        OPF_ALG_DC = 500;
    case 'MOSEK'
        OPF_ALG_DC = 600;
    case 'GUROBI'
        OPF_ALG_DC = 700;
end

%% KNITRO_OPT
if isfield(opt_s, 'knitro')
    KNITRO_OPT  = opt_s.knitro.opt;
else
    KNITRO_OPT  = old(58);
end

%% IPOPT_OPT
if isfield(opt_s, 'ipopt')
    IPOPT_OPT  = opt_s.ipopt.opt;
else
    IPOPT_OPT  = old(58);
end

%% MINOPF options
if isfield(opt_s, 'minopf')
    MINOPF_OPTS = [
        opt_s.minopf.feastol;   %% 61 - MNS_FEASTOL
        opt_s.minopf.rowtol;    %% 62 - MNS_ROWTOL
        opt_s.minopf.xtol;      %% 63 - MNS_XTOL
        opt_s.minopf.majdamp;   %% 64 - MNS_MAJDAMP
        opt_s.minopf.mindamp;   %% 65 - MNS_MINDAMP
        opt_s.minopf.penalty;   %% 66 - MNS_PENALTY_PARM
        opt_s.minopf.major_it;  %% 67 - MNS_MAJOR_IT
        opt_s.minopf.minor_it;  %% 68 - MNS_MINOR_IT
        opt_s.minopf.max_it;    %% 69 - MNS_MAX_IT
        opt_s.minopf.verbosity; %% 70 - MNS_VERBOSITY
        opt_s.minopf.core;      %% 71 - MNS_CORE
        opt_s.minopf.supbasic_lim;  %% 72 - MNS_SUPBASIC_LIM
        opt_s.minopf.mult_price;%% 73 - MNS_MULT_PRICE
    ];
else
    MINOPF_OPTS = old(61:73);
end

%% FORCE_PC_EQ_P0
if isfield(opt_s, 'sopf') && isfield(opt_s.sopf, 'force_Pc_eq_P0')
    FORCE_PC_EQ_P0 = opt_s.sopf.force_Pc_eq_P0;
else
    FORCE_PC_EQ_P0 = 0;
end

%% PDIPM options
if isfield(opt_s, 'pdipm')
    PDIPM_OPTS = [
        opt_s.pdipm.feastol;    %% 81 - PDIPM_FEASTOL
        opt_s.pdipm.gradtol;    %% 82 - PDIPM_GRADTOL
        opt_s.pdipm.comptol;    %% 83 - PDIPM_COMPTOL
        opt_s.pdipm.costtol;    %% 84 - PDIPM_COSTTOL
        opt_s.pdipm.max_it;     %% 85 - PDIPM_MAX_IT
        opt_s.pdipm.sc.red_it;  %% 86 - SCPDIPM_RED_IT
    ];
else
    PDIPM_OPTS = old(81:86);
end

%% TRALM options
if isfield(opt_s, 'tralm')
    TRALM_OPTS = [
        opt_s.tralm.feastol;    %% 87 - TRALM_FEASTOL
        opt_s.tralm.primaltol;  %% 88 - TRALM_PRIMETOL
        opt_s.tralm.dualtol;    %% 89 - TRALM_DUALTOL
        opt_s.tralm.costtol;    %% 90 - TRALM_COSTTOL
        opt_s.tralm.major_it;   %% 91 - TRALM_MAJOR_IT
        opt_s.tralm.minor_it;   %% 92 - TRALM_MINOR_IT
    ];
else
    TRALM_OPTS = old(87:92);
end

%% SMOOTHING_RATIO
if strcmp(upper(opt_s.opf.ac.solver), 'TRALM') && isfield(opt_s, 'tralm')
    SMOOTHING_RATIO = opt_s.tralm.smooth_ratio;
elseif isfield(opt_s, 'pdipm')
    SMOOTHING_RATIO = opt_s.pdipm.sc.smooth_ratio;
else
    SMOOTHING_RATIO = old(93);
end

%% CPLEX options
if isfield(opt_s, 'cplex')
    CPLEX_OPTS = [
        opt_s.cplex.lpmethod;   %% 95 - CPLEX_LPMETHOD
        opt_s.cplex.qpmethod;   %% 96 - CPLEX_QPMETHOD
        opt_s.cplex.opt;        %% 97 - CPLEX_OPT
    ];
else
    CPLEX_OPTS = old(95:97);
end

%% MOSEK options
if isfield(opt_s, 'mosek')
    MOSEK_OPTS = [
        opt_s.mosek.lp_alg;     %% 111 - MOSEK_LP_ALG
        opt_s.mosek.max_it;     %% 112 - MOSEK_MAX_IT
        opt_s.mosek.gap_tol;    %% 113 - MOSEK_GAP_TOL
        opt_s.mosek.max_time;   %% 114 - MOSEK_MAX_TIME
        opt_s.mosek.num_threads;%% 115 - MOSEK_NUM_THREADS
        opt_s.mosek.opt;        %% 116 - MOSEK_OPT
    ];
else
    MOSEK_OPTS = old(111:116);
end

%% Gurobi options
if isfield(opt_s, 'gurobi')
    GUROBI_OPTS = [
        opt_s.gurobi.method;    %% 121 - GRB_METHOD
        opt_s.gurobi.timelimit; %% 122 - GRB_TIMELIMIT
        opt_s.gurobi.threads;   %% 123 - GRB_THREADS
        opt_s.gurobi.opt;       %% 124 - GRB_OPT
    ];
else
    GUROBI_OPTS = old(121:124);
end

opt_v = [
        %% power flow options
        PF_ALG;                 %% 1  - PF_ALG
        opt_s.pf.tol;           %% 2  - PF_TOL
        opt_s.pf.nr.max_it;     %% 3  - PF_MAX_IT
        opt_s.pf.fd.max_it;     %% 4  - PF_MAX_IT_FD
        opt_s.pf.gs.max_it;     %% 5  - PF_MAX_IT_GS
        opt_s.pf.enforce_q_lims;%% 6  - ENFORCE_Q_LIMS
        0;                      %% 7  - RESERVED7
        0;                      %% 8  - RESERVED8
        0;                      %% 9  - RESERVED9
        PF_DC;                  %% 10 - PF_DC
        
        %% OPF options
        OPF_ALG;                %% 11 - OPF_ALG
        0;                      %% 12 - RESERVED12 (was OPF_ALG_POLY = 100)
        0;                      %% 13 - RESERVED13 (was OPF_ALG_PWL = 200)
        0;                      %% 14 - RESERVED14 (was OPF_POLY2PWL_PTS = 10)
        0;                      %% 15 - OPF_NEQ (removed)
        opt_s.opf.violation;    %% 16 - OPF_VIOLATION
        CONSTR_TOL_X;           %% 17 - CONSTR_TOL_X
        CONSTR_TOL_F;           %% 18 - CONSTR_TOL_F
        CONSTR_MAX_IT;          %% 19 - CONSTR_MAX_IT
        old(20);                %% 20 - LPC_TOL_GRAD (removed)
        old(21);                %% 21 - LPC_TOL_X (removed)
        old(22);                %% 22 - LPC_MAX_IT (removed)
        old(23);                %% 23 - LPC_MAX_RESTART (removed)
        OPF_FLOW_LIM;           %% 24 - OPF_FLOW_LIM
        opt_s.opf.ignore_angle_lim; %% 25 - OPF_IGNORE_ANG_LIM
        OPF_ALG_DC;             %% 26 - OPF_ALG_DC
        0;                      %% 27 - RESERVED27
        0;                      %% 28 - RESERVED28
        0;                      %% 29 - RESERVED29
        0;                      %% 30 - RESERVED30
        
        %% output options
        opt_s.verbose;          %% 31 - VERBOSE
        opt_s.out.all;          %% 32 - OUT_ALL
        opt_s.out.sys_sum;      %% 33 - OUT_SYS_SUM
        opt_s.out.area_sum;     %% 34 - OUT_AREA_SUM
        opt_s.out.bus;          %% 35 - OUT_BUS
        opt_s.out.branch;       %% 36 - OUT_BRANCH
        opt_s.out.gen;          %% 37 - OUT_GEN
        opt_s.out.lim.all;      %% 38 - OUT_ALL_LIM
        opt_s.out.lim.v;        %% 39 - OUT_V_LIM
        opt_s.out.lim.line;     %% 40 - OUT_LINE_LIM
        opt_s.out.lim.pg;       %% 41 - OUT_PG_LIM
        opt_s.out.lim.qg;       %% 42 - OUT_QG_LIM
        0;                      %% 43 - RESERVED43 (was OUT_RAW)
        opt_s.out.force;        %% 44 - OUT_FORCE
        0;                      %% 45 - RESERVED45
        0;                      %% 46 - RESERVED46
        0;                      %% 47 - RESERVED47
        0;                      %% 48 - RESERVED48
        0;                      %% 49 - RESERVED49
        0;                      %% 50 - RESERVED50
        
        %% other options
        old(51);                %% 51 - SPARSE_QP (removed)
        opt_s.opf.return_raw_der;   %% 52 - RETURN_RAW_DER
        0;                      %% 53 - RESERVED53
        0;                      %% 54 - RESERVED54
        FMC_ALG;                %% 55 - FMC_ALG
        0;                      %% 56 - RESERVED56
        0;                      %% 57 - RESERVED57
        KNITRO_OPT;             %% 58 - KNITRO_OPT
        0;                      %% 59 - RESERVED59
        IPOPT_OPT;              %% 60 - IPOPT_OPT
        
        %% MINOPF options
        MINOPF_OPTS;            %% 61-73 - MNS_FEASTOL-MNS_MULT_PRICE
        0;                      %% 74 - RESERVED74
        0;                      %% 75 - RESERVED75
        0;                      %% 76 - RESERVED76
        0;                      %% 77 - RESERVED77
        0;                      %% 78 - RESERVED78
        0;                      %% 79 - RESERVED79
        FORCE_PC_EQ_P0;         %% 80 - FORCE_PC_EQ_P0, for c3sopf
        
        %% MIPS, PDIPM, SC-PDIPM, and TRALM options
        PDIPM_OPTS;             %% 81-86 - PDIPM_FEASTOL-SCPDIPM_RED_IT
        TRALM_OPTS;             %% 87-92 - TRALM_FEASTOL-TRALM_MINOR_IT
        SMOOTHING_RATIO;        %% 93 - SMOOTHING_RATIO
        0;                      %% 94 - RESERVED94
        
        %% CPLEX options
        CPLEX_OPTS;             %% 95-97 - CPLEX_LPMETHOD-CPLEX_OPT
        0;                      %% 98 - RESERVED98
        0;                      %% 99 - RESERVED99
        0;                      %% 100 - RESERVED100
        0;                      %% 101 - RESERVED101
        0;                      %% 102 - RESERVED102
        0;                      %% 103 - RESERVED103
        0;                      %% 104 - RESERVED104
        0;                      %% 105 - RESERVED105
        0;                      %% 106 - RESERVED106
        0;                      %% 107 - RESERVED107
        0;                      %% 108 - RESERVED108
        0;                      %% 109 - RESERVED109
        0;                      %% 110 - RESERVED110

        %% MOSEK options
        MOSEK_OPTS;             %% 111-116 - MOSEK_LP_ALG-MOSEK_OPT
        0;                      %% 117 - RESERVED117
        0;                      %% 118 - RESERVED118
        0;                      %% 119 - RESERVED119
        0;                      %% 120 - RESERVED120

        %% Gurobi options
        GUROBI_OPTS;            %% 121-124 - GRB_METHOD-GRB_OPT
    ];


%%-------------------------------------------------------------------
function optt = mpoption_default()
if DEBUG, fprintf('mpoption_default()\n'); end
persistent opt;             %% cache this for speed
if ~isstruct(opt)
    opt = struct(...
        'v',                    mpoption_version, ...   %% version
        'model',                'AC', ...
        'pf',                   struct(...
            'alg',                  'NR', ...
            'tol',                  1e-8, ...
            'nr',                   struct(...
                'max_it',               10  ), ...
            'fd',                   struct(...
                'max_it',               30  ), ...
            'gs',                   struct(...
                'max_it',               1000  ), ...
            'enforce_q_lims',       0   ), ...
        'cpf',                  struct(...
            'parameterization',     3, ...
            'stop_at',              'NOSE', ...     %% 'NOSE', <lam val>, 'FULL'
            'enforce_p_lims',       0, ...
            'enforce_q_lims',       0, ...
            'step',                 0.05, ...
            'step_min',             1e-4, ...
            'step_max',             0.2, ...
            'adapt_step',           0, ...
            'adapt_step_damping',   0.7, ...
            'adapt_step_tol',       1e-3, ...
            'target_lam_tol',       1e-5, ...
            'nose_tol',             1e-5, ...
            'p_lims_tol',           0.01, ...
            'q_lims_tol',           0.01, ...
            'plot',                 struct(...
                'level',                0, ...
                'bus',                  []  ), ...
            'user_callback',        ''  ), ...
        'opf',                  struct(...
            'ac',                   struct(...
                'solver',               'DEFAULT'   ), ...
            'dc',                   struct(...
                'solver',               'DEFAULT'   ), ...
            'violation',            5e-6, ...
            'use_vg',               0, ...
            'flow_lim',             'S', ...
            'ignore_angle_lim',     0, ...
            'init_from_mpc',        -1, ...
            'return_raw_der',       0   ), ...
        'verbose',              1, ...
        'out',                  struct(...
            'all',                  -1, ...
            'sys_sum',              1, ...
            'area_sum',             0, ...
            'bus',                  1, ...
            'branch',               1, ...
            'gen',                  0, ...
            'lim',                  struct(...
                'all',                  -1, ...
                'v',                    1, ...
                'line',                 1, ...
                'pg',                   1, ...
                'qg',                   1   ), ...
            'force',                0, ...
            'suppress_detail',      -1  ), ...
        'mips',                 struct(...  %% see mpoption_info_mips() for optional fields
            'step_control',         0, ...
            'linsolver',            '', ...
            'feastol',              0, ...
            'gradtol',              1e-6, ...
            'comptol',              1e-6, ...
            'costtol',              1e-6, ...
            'max_it',               150, ...
            'sc',                   struct(...
                'red_it',               20  )), ...
        'exp',                  struct(... %% experimental options
            'sys_wide_zip_loads',   struct(...
                'pw',                   [], ...
                'qw',                   []  )) ...
    );
    opt_pkgs = mpoption_optional_pkgs();
    for k = 1:length(opt_pkgs)
        fname = ['mpoption_info_' opt_pkgs{k}];
        if exist(fname, 'file') == 2
            opt = nested_struct_copy(opt, feval(fname, 'D'));
        end
    end
end
optt = opt;

%%-------------------------------------------------------------------
function optt = mpoption_optional_fields()
if DEBUG, fprintf('mpoption_optional_fields()\n'); end
persistent opt;         %% cache this for speed
if ~isstruct(opt)
    opt_pkgs = mpoption_optional_pkgs();
    opt = struct;
    for k = 1:length(opt_pkgs)
        fname = ['mpoption_info_' opt_pkgs{k}];
        if exist(fname, 'file') == 2
            opt = nested_struct_copy(opt, feval(fname, 'V'));
        end
    end
end
optt = opt;

%% globals
%%-------------------------------------------------------------------
function v = mpoption_version
v = 12;     %% version number of MATPOWER options struct
            %% (must be incremented every time structure is updated)
            %% v1   - first version based on struct (MATPOWER 5.0b1)
            %% v2   - added 'linprog' and 'quadprog' fields
            %% v3   - (forgot to increment v) added 'out.suppress_detail'
            %%        field
            %% v4   - (forgot to increment v) MIPS 1.1, added optional
            %%        fields to 'mips' options: xi, sigma, z0, alpha_min,
            %%        rho_min, rho_max, mu_threshold and max_stepsize
            %% v5   - (forgot to increment v) added 'opf.init_from_mpc'
            %%        field (MATPOWER 5.0)
            %% v6   - added 'clp' field
            %% v7   - added 'intlinprog' field
            %% v8   - MIPS 1.2, added 'linsolver' field to
            %%        'mips' options
            %% v9   - added 'exp' for experimental fields, specifically
            %%        'sys_wide_zip_loads.pw', 'sys_wide_zip_loads.qw'
            %% v10  - added 'most' field
            %% v11  - added 'cpf' options 'adapt_step_damping',
            %%        'enforce_p_lims', 'enforce_q_lims', 'target_lam_tol'
            %%        'nose_tol', 'p_lims_tol' and 'q_lims_tol',
            %%        removed option 'cpf.user_callback_args'
            %% v12  - added option 'opf.use_vg'

%%-------------------------------------------------------------------
function db_level = DEBUG
db_level = 0;

%%-------------------------------------------------------------------
function pkgs = mpoption_optional_pkgs()
pkgs = {...
    'clp', 'cplex', 'fmincon', 'gurobi', 'glpk', 'intlinprog', 'ipopt', ...
    'knitro', 'linprog', 'minopf', 'most', 'mosek', 'quadprog', 'sdp_pf', ...
    'sopf', 'tspopf', 'yalmip' ...
};


function [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas)
%EXT2INT   Converts external to internal indexing.
%
%   This function has two forms, (1) the old form that operates on
%   and returns individual matrices and (2) the new form that operates
%   on and returns an entire MATPOWER case struct.
%
%   1.  [I2E, BUS, GEN, BRANCH, AREAS] = EXT2INT(BUS, GEN, BRANCH, AREAS)
%       [I2E, BUS, GEN, BRANCH] = EXT2INT(BUS, GEN, BRANCH)
%
%   If the first argument is a matrix, it simply converts from (possibly
%   non-consecutive) external bus numbers to consecutive internal bus
%   numbers which start at 1. Changes are made to BUS, GEN and BRANCH,
%   which are returned along with a vector of indices I2E that can be
%   passed to INT2EXT to perform the reverse conversion, where
%   EXTERNAL_BUS_NUMBER = I2E(INTERNAL_BUS_NUMBER).
%   AREAS is completely ignored and is only included here for backward
%   compatibility of the API.
%
%   Examples:
%       [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);
%       [i2e, bus, gen, branch] = ext2int(bus, gen, branch);
%
%   2.  MPC = EXT2INT(MPC)
%
%   If the input is a single MATPOWER case struct, then all isolated
%   buses, off-line generators and branches are removed along with any
%   generators or branches connected to isolated buses. Then the
%   buses are renumbered consecutively, beginning at 1, and the
%   generators are sorted by increasing bus number. Any 'ext2int'
%   callback routines registered in the case are also invoked
%   automatically. All of the related indexing information and the
%   original data matrices are stored in an 'order' field in the struct
%   to be used by INT2EXT to perform the reverse conversions. If the
%   case is already using internal numbering it is returned unchanged.
%
%   Example:
%       mpc = ext2int(mpc);
%
%   The 'order' field of MPC used to store the indexing information
%   needed for subsequent internal to external conversion is structured
%   as:
%
%       order
%           state       'i' | 'e'
%           ext | int
%               bus
%               branch
%               gen
%               gencost
%               A
%               N
%           bus
%               e2i
%               i2e
%               status
%                   on
%                   off
%           gen
%               e2i
%               i2e
%               status
%                   on
%                   off
%           branch
%               status
%                   on
%                   off
%
%   See also INT2EXT, E2I_FIELD, E2I_DATA.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isstruct(bus)
    mpc = bus;
    if nargin == 1
        first = ~isfield(mpc, 'order');
        if first || mpc.order.state == 'e'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS] = idx_brch;

            %% initialize order
            if first
                status = struct('on',   [], ...
                                'off',  []  );
                tmp = struct( ...
                        'e2i',      [], ...
                        'i2e',      [], ...
                        'status',   status ...
                    );
                o = struct( ...
                        'ext',      struct( ...
                                'bus',      [], ...
                                'branch',   [], ...
                                'gen',      [] ...
                            ), ...
                        'bus',      tmp, ...
                        'gen',      tmp, ...
                        'branch',   struct('status', status) ...
                    );
            else
                o = mpc.order;
            end

            %% sizes
            nb = size(mpc.bus, 1);
            ng = size(mpc.gen, 1);
            ng0 = ng;
            if isfield(mpc, 'A') && size(mpc.A, 2) < 2*nb + 2*ng
                dc = 1;
            elseif isfield(mpc, 'N') && size(mpc.N, 2) < 2*nb + 2*ng
                dc = 1;
            else
                dc = 0;
            end

            %% save data matrices with external ordering
            o.ext.bus    = mpc.bus;
            o.ext.branch = mpc.branch;
            o.ext.gen    = mpc.gen;

            %% check that all buses have a valid BUS_TYPE
            bt = mpc.bus(:, BUS_TYPE);
            err = find(~(bt == PQ | bt == PV | bt == REF | bt == NONE));
            if ~isempty(err)
                error('ext2int: bus %d has an invalid BUS_TYPE', err);
            end

            %% determine which buses, branches, gens are connected & in-service
            n2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
            bs = (bt ~= NONE);                      %% bus status
            o.bus.status.on     = find(  bs );      %% connected
            o.bus.status.off    = find( ~bs );      %% isolated
            gs = ( mpc.gen(:, GEN_STATUS) > 0 & ... %% gen status
                    bs(n2i(mpc.gen(:, GEN_BUS))) );
            o.gen.status.on     = find(  gs );      %% on and connected
            o.gen.status.off    = find( ~gs );      %% off or isolated
            brs = ( mpc.branch(:, BR_STATUS) & ...  %% branch status
                    bs(n2i(mpc.branch(:, F_BUS))) & ...
                    bs(n2i(mpc.branch(:, T_BUS))) );
            o.branch.status.on  = find(  brs );     %% on and connected
            o.branch.status.off = find( ~brs );

            %% delete stuff that is "out"
            if ~isempty(o.bus.status.off)
                mpc.bus(o.bus.status.off, :) = [];
            end
            if ~isempty(o.branch.status.off)
                mpc.branch(o.branch.status.off, :) = [];
            end
            if ~isempty(o.gen.status.off)
                mpc.gen(o.gen.status.off, :) = [];
            end

            %% update size
            nb = size(mpc.bus, 1);

            %% apply consecutive bus numbering
            o.bus.i2e = mpc.bus(:, BUS_I);
            o.bus.e2i = sparse(max(o.bus.i2e), 1);
            o.bus.e2i(o.bus.i2e) = (1:nb)';
            mpc.bus(:, BUS_I)       = o.bus.e2i( mpc.bus(:, BUS_I)      );
            mpc.gen(:, GEN_BUS)     = o.bus.e2i( mpc.gen(:, GEN_BUS)    );
            mpc.branch(:, F_BUS)    = o.bus.e2i( mpc.branch(:, F_BUS)   );
            mpc.branch(:, T_BUS)    = o.bus.e2i( mpc.branch(:, T_BUS)   );

            %% reorder gens in order of increasing bus number
            [tmp, o.gen.e2i] = sort(mpc.gen(:, GEN_BUS));
            [tmp, o.gen.i2e] = sort(o.gen.e2i);
            mpc.gen = mpc.gen(o.gen.e2i, :);

            if isfield(o, 'int')
                o = rmfield(o, 'int');
            end
            o.state = 'i';
            mpc.order = o;

            %% update gencost, A and N
            if isfield(mpc, 'gencost')
                ordering = {'gen'};         %% Pg cost only
                if size(mpc.gencost, 1) == 2*ng0
                    ordering{2} = 'gen';    %% include Qg cost
                end
                mpc = e2i_field(mpc, 'gencost', ordering);
            end
            if isfield(mpc, 'A') || isfield(mpc, 'N')
                if dc
                    ordering = {'bus', 'gen'};
                else
                    ordering = {'bus', 'bus', 'gen', 'gen'};
                end
            end
            if isfield(mpc, 'A')
                mpc = e2i_field(mpc, 'A', ordering, 2);
            end
            if isfield(mpc, 'N')
                mpc = e2i_field(mpc, 'N', ordering, 2);
            end

            %% execute userfcn callbacks for 'ext2int' stage
            if isfield(mpc, 'userfcn')
                mpc = run_userfcn(mpc.userfcn, 'ext2int', mpc);
            end
        end

        i2e = mpc;
    else                    %% convert extra data
        ordering = branch;              %% rename argument
        if nargin < 4
            dim = 1;
        else
            dim = areas;                %% rename argument
        end
        if ischar(gen) || iscell(gen)   %% field
            warning('Calls of the form MPC = EXT2INT(MPC, ''FIELD_NAME'', ...) have been deprecated. Please replace EXT2INT with E2I_FIELD.');
            i2e = e2i_field(mpc, gen, branch, dim);
        else                            %% value
            warning('Calls of the form VAL = EXT2INT(MPC, VAL, ...) have been deprecated. Please replace EXT2INT with E2I_DATA.');
            i2e = e2i_data(mpc, gen, branch, dim);
        end
    end
else            %% old form
    %% define names for columns to data matrices
    [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    [GEN_BUS] = idx_gen;
    [F_BUS, T_BUS] = idx_brch;

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';

    %% renumber buses consecutively
    bus(:, BUS_I)               = e2i( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = e2i( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = e2i( branch(:, F_BUS)         );
    branch(:, T_BUS)            = e2i( branch(:, T_BUS)         );
end

function mpc = e2i_field(mpc, field, ordering, dim)
%E2I_FIELD   Converts fields of MPC from external to internal indexing.
%
%   This function performs several different tasks, depending on the
%   arguments passed.
%
%   MPC = E2I_FIELD(MPC, FIELD, ORDERING)
%   MPC = E2I_FIELD(MPC, FIELD, ORDERING, DIM)
%
%   When given a case struct that has already been converted to
%   internal indexing, this function can be used to convert other data
%   structures as well by passing in 2 or 3 extra parameters in
%   addition to the case struct.
%
%   The 2nd argument is a string or cell array of strings, specifying
%   a field in the case struct whose value should be converted by
%   a corresponding call to E2I_DATA. The field can contain either a
%   numeric or a cell array. The converted value is stored back in the
%   specified field, the original value is saved for later use and the
%   updated case struct is returned. If FIELD is a cell array of strings,
%   they specify nested fields.
%
%   The 3rd and optional 4th arguments are simply passed along to
%   the call to E2I_DATA.
%
%   Examples:
%       mpc = e2i_field(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match internal generator
%       ordering.
%
%       mpc = e2i_field(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match internal
%       generator ordering.
%
%   See also I2E_FIELD, E2I_DATA, EXT2INT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    dim = 1;
end
if ischar(field)
    mpc.order.ext.(field) = mpc.(field);
    mpc.(field) = e2i_data(mpc, mpc.(field), ordering, dim);
else    %% iscell(field)
    for k = 1:length(field)
        s(k).type = '.';
        s(k).subs = field{k};
    end
    mpc.order.ext = subsasgn(mpc.order.ext, s, subsref(mpc, s));
    mpc = subsasgn(mpc, s, ...
        e2i_data(mpc, subsref(mpc, s), ordering, dim) );
end

function newval = e2i_data(mpc, val, ordering, dim)
%E2I_DATA   Converts data from external to internal indexing.
%
%   VAL = E2I_DATA(MPC, VAL, ORDERING)
%   VAL = E2I_DATA(MPC, VAL, ORDERING, DIM)
%
%   When given a case struct that has already been converted to
%   internal indexing, this function can be used to convert other data
%   structures as well by passing in 2 or 3 extra parameters in
%   addition to the case struct. If the value passed in the 2nd
%   argument is a column vector or cell array, it will be converted
%   according to the ORDERING specified by the 3rd argument (described
%   below). If VAL is an n-dimensional matrix or cell array, then the
%   optional 4th argument (DIM, default = 1) can be used to specify
%   which dimension to reorder. The return value in this case is the
%   value passed in, converted to internal indexing.
%
%   The 3rd argument, ORDERING, is used to indicate whether the data
%   corresponds to bus-, gen- or branch-ordered data. It can be one
%   of the following three strings: 'bus', 'gen' or 'branch'. For
%   data structures with multiple blocks of data, ordered by bus,
%   gen or branch, they can be converted with a single call by
%   specifying ORDERING as a cell array of strings.
%
%   Any extra elements, rows, columns, etc. beyond those indicated
%   in ORDERING, are not disturbed.
%
%   Examples:
%       A_int = e2i_data(mpc, A_ext, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       external to internal ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_int = e2i_data(mpc, gencost_ext, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%   See also I2E_DATA, E2I_FIELD, EXT2INT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if ~isfield(mpc, 'order')
    error('e2i_data: mpc does not have the ''order'' field required to convert from external to internal numbering.');
end
o = mpc.order;
if o.state ~= 'i'
    error('e2i_data: mpc does not have internal ordering data available, call ext2int first');
end
if nargin < 4
    dim = 1;
end
if ischar(ordering)         %% single set
    if strcmp(ordering, 'gen')
        idx = o.(ordering).status.on(o.(ordering).e2i);
    else
        idx = o.(ordering).status.on;
    end
    newval = get_reorder(val, idx, dim);
else                        %% multiple sets
    b = 0;  %% base
    for k = 1:length(ordering)
        n = size(o.ext.(ordering{k}), 1);
        v = get_reorder(val, b+(1:n), dim);
        new_v{k} = e2i_data(mpc, v, ordering{k}, dim);
        b = b + n;
    end
    n = size(val, dim);
    if n > b                %% the rest
        v = get_reorder(val, b+1:n, dim);
        new_v{length(new_v)+1} = v;
    end
    newval = cat(dim, new_v{:});
end

function B = get_reorder(A, idx, dim)
%GET_REORDER    Returns A with one of its dimensions indexed.
%
%   B = GET_REORDER(A, IDX, DIM)
%
%   Returns A(:, ..., :, IDX, :, ..., :), where DIM determines
%   in which dimension to place the IDX.
%
%   See also SET_REORDER.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

ndim = ndims(A);
s.type = '()';
s.subs = cell(1, ndim);
for k = 1:ndim
    if k == dim
        s.subs{k} = idx;
    else
        s.subs{k} = ':';
    end
end
B = subsref(A, s);

function [loss, fchg, tchg, dloss_dV, dchg_dVm] = get_losses(baseMVA, bus, branch)
%GET_LOSSES   Returns series losses (and reactive injections) per branch.
%
%   LOSS = GET_LOSSES(RESULTS)
%   LOSS = GET_LOSSES(BASEMVA, BUS, BRANCH)
%
%   [LOSS, CHG] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG, DLOSS_DV] = GET_LOSSES(RESULTS)
%   [LOSS, FCHG, TCHG, DLOSS_DV, DCHG_DVM] = GET_LOSSES(RESULTS)
%
%   Computes branch series losses, and optionally reactive injections from
%   line charging, as functions of bus voltages and branch parameters, using the
%   following formulae:
%
%       loss = abs( Vf / tau - Vt ) ^ 2 / (Rs - j Xs)
%       fchg = abs( Vf / tau ) ^ 2 * Bc / 2
%       tchg = abs( Vt ) ^ 2 * Bc / 2
%
%   Optionally, computes the partial derivatives of the line losses with
%   respect to voltage angles and magnitudes.
%
%   Input:
%       RESULTS - a MATPOWER case struct with bus voltages corresponding to
%                 a valid power flow solution.
%                 (Can optionally be specified as individual fields BASEMVA,
%                  BUS, and BRANCH.)
%
%   Output(s):
%       LOSS - complex NL x 1 vector of losses (in MW), where NL is the number
%              of branches in the system, representing only the losses in the
%              series impedance element of the PI model for each branch.
%       CHG -  NL x 1 vector of total reactive injection for each line
%              (in MVAr), representing the line charging injections of both
%              of the shunt elements of PI model for each branch.
%       FCHG - Same as CHG, but for the element at the "from" end of the
%              branch only.
%       TCHG - Same as CHG, but for the element at the "to" end of the branch.
%       DLOSS_DV - Struct with partial derivatives of LOSS with respect to bus
%              voltages, with fields:
%           .a  - Partial with respect to bus voltage angles.
%           .m  - Partial with respect to bus voltage magnitudes.
%       DCHG_DVM - Struct with partial derivatives of FCHG and TCHG with
%              respect to bus voltage magnitudes, with fields:
%           .f  - Partial of FCHG with respect to bus voltage magnitudes.
%           .t  - Partial of TCHG with respect to bus voltage magnitudes.
%
%   Example:
%       results = runpf(mycase);
%       total_system_real_losses = sum(real(get_losses(results)));
%
%       [loss, fchg, tchg, dloss_dV] = get_losses(results);

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if isstruct(baseMVA)
    mpc = baseMVA;
    [baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% create map of external bus numbers to bus indices
i2e = bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(bus, 1))';
out = find(branch(:, BR_STATUS) == 0);          %% out-of-service branches

%% sizes of things
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches

%% construct complex bus voltage vector
V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));

%% parameters
Cf = sparse(1:nl, e2i(branch(:, F_BUS)), branch(:, BR_STATUS), nl, nb);
Ct = sparse(1:nl, e2i(branch(:, T_BUS)), branch(:, BR_STATUS), nl, nb);
tap = ones(nl, 1);                              %% default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));                    %% indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  %% include transformer tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
A = spdiags(1 ./ tap, 0, nl, nl) * Cf - Ct;
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X));
Vdrop = A * V;      %% vector of voltage drop across series impedance element
loss = baseMVA * Ysc .* Vdrop .* conj(Vdrop);
% loss = baseMVA * abs(V(e2i(branch(:, F_BUS))) ./ tap - V(e2i(branch(:, T_BUS)))) .^ 2 ./ ...
%             (branch(:, BR_R) - 1j * branch(:, BR_X));
% loss(out) = 0;

if nargout > 1
    Vf = Cf * V;
    Vt = Ct * V;
    fchg = real(baseMVA / 2 * branch(:, BR_B) .* Vf .* conj(Vf) ./ (tap .* conj(tap)));
    tchg = real(baseMVA / 2 * branch(:, BR_B) .* Vt .* conj(Vt));
%     fchg = abs(V(e2i(branch(:, F_BUS))) ./ tap) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
%     tchg = abs(V(e2i(branch(:, T_BUS)))       ) .^ 2 .* branch(:, BR_B) * baseMVA / 2;
    fchg(out) = 0;
    tchg(out) = 0;

    if nargout == 2
        fchg = fchg + tchg;
    end
end

if nargout > 3
    B = spdiags(A * V, 0, nl, nl) * conj(A) * spdiags(conj(V), 0, nb, nb);
    dYsc = spdiags(Ysc, 0, nl, nl);
    dloss_dV = struct(...
        'a', -1j * baseMVA * dYsc * (B - conj(B)), ...
        'm',       baseMVA * dYsc * (B + conj(B)) * spdiags(1 ./ abs(V), 0, nb, nb) ...
    );
    if nargout > 4
        Bc = spdiags(branch(:, BR_B), 0, nl, nl);
        tt = spdiags(1 ./ (tap .* conj(tap)), 0, nl, nl);
        dchg_dVm = struct(...
            'f', baseMVA * Bc * tt * spdiags(Cf * bus(:, VM), 0, nb, nb) * Cf, ...
            't', baseMVA * Bc      * spdiags(Ct * bus(:, VM), 0, nb, nb) * Ct ...
        );
    end
end

function rv = mpver(varargin)
%MPVER  Prints or returns MATPOWER version info for current installation.
%   V = MPVER returns the current MATPOWER version number.
%   V = MPVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MPVER without assigning the
%   return value prints the version and release date of the current
%   installation of MATPOWER, MATLAB (or Octave), the Optimization Toolbox,
%   MIPS and any optional MATPOWER packages.

%   MATPOWER
%   Copyright (c) 2005-2015, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% the following only works if MATPOWER is explicitly in the path,
%% but not if it is only in the current working directory
% fs = filesep;
% p = fileparts(which('runpf'));
% if ~strcmp(p(1),fs)
%   [t, p] = strtok(p, filesep);
% end
% p = p(2:end);
% v{1} = ver(p);

v{1} = struct(  'Name',     'MATPOWER', ... 
                'Version',  '6.0', ...
                'Release',  '', ...
                'Date',     '16-Dec-2016' );
if nargout > 0
    if nargin > 0
        rv = v{1};
    else
        rv = v{1}.Version;
    end
else
    if have_fcn('octave')
        v{2} = ver('octave');
    else
        v{2} = ver('matlab');
        if length(v{2}) > 1
            warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
            v{2} = v{2}(1);
        end
    end
    v{3} = ver('optim');
    if length(v{3}) > 1
        warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''optim'' on your path. Check each element of the output of ver(''optim'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
        v{3} = v{3}(1);
    end
    for n = 1:3
        if n == 3
            if isempty(v{3})
                fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
                continue;
            elseif have_fcn('matlab') && ~license('test', 'optimization_toolbox')
                fprintf('\n%-22s -- no license --', 'Optimization Toolbox');
                continue;
            end
        end
        fprintf('\n%-22s Version %-9s', v{n}.Name, v{n}.Version);
        if ~isempty(v{n}.Date)
            fprintf('  %11s', v{n}.Date);
            if ~isempty(v{n}.Release)
                fprintf('   Release: %-10s', v{n}.Release);
            end
        end
    end
    fprintf('\n');
    mipsver;
    if have_fcn('most')
        mostver;
    else
        fprintf('%-22s -- not installed --\n', 'MOST');
    end
    if have_fcn('sdp_pf')
        sdp_pf_ver;
    else
        fprintf('%-22s -- not installed --\n', 'SDP_PF');
    end
    if have_fcn('yalmip')
        s = have_fcn('yalmip', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'YALMIP', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'YALMIP');
    end
    if have_fcn('bpmpd')
        if exist('bpver', 'file') == 2
            bpver;
        else
            fprintf('%-22s Version 2.21 or earlier\n', 'BPMPD_MEX');
        end
    else
        fprintf('%-22s -- not installed --\n', 'BPMPD_MEX');
    end
    if have_fcn('clp')
        s = have_fcn('clp', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'CLP', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CLP');
    end
    if have_fcn('cplex')
        s = have_fcn('cplex', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'CPLEX', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CPLEX');
    end
    if have_fcn('glpk')
        s = have_fcn('glpk', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'GLPK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'GLPK');
    end
    gurobiver;
    if have_fcn('ipopt')
        s = have_fcn('ipopt', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'IPOPT', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'IPOPT');
    end
    if have_fcn('knitro')
        s = have_fcn('knitro', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'KNITRO', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'KNITRO');
    end
    if have_fcn('minopf')
        if exist('minopfver', 'file') == 2
            minopfver;
        else
            fprintf('%-22s Version 3.0b2 or earlier\n', 'MINOPF');
        end
    else
        fprintf('%-22s -- not installed --\n', 'MINOPF');
    end
    if have_fcn('mosek')
        s = have_fcn('mosek', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'MOSEK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'MOSEK');
    end
    if have_fcn('pardiso')
        s = have_fcn('pardiso', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'PARDISO', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'PARDISO');
    end
    if have_fcn('pdipmopf')
        pdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'PDIPMOPF');
    end
    if have_fcn('scpdipmopf')
        scpdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'SCPDIPMOPF');
    end
    if have_fcn('sdpt3')
        s = have_fcn('sdpt3', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SDPT3', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SDPT3');
    end
    if have_fcn('sedumi')
        s = have_fcn('sedumi', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SeDuMi', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SeDuMi');
    end
    if have_fcn('tralmopf')
        tralmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'TRALMOPF');
    end

    fprintf('%-22s %s\n\n', 'Architecture:', computer);
    
    fprintf('  MATPOWER %s is distributed under the 3-clause BSD License.\n', v{1}.Version);
    fprintf('  Please see the LICENSE file for details.\n\n');
end

function Sd = makeSdzip(baseMVA, bus, mpopt)
%MAKESDZIP   Builds vectors of nominal complex bus power demands for ZIP loads.
%   SD = MAKESDZIP(BASEMVA, BUS, MPOPT) returns a struct with three fields,
%   each an nb x 1 vectors. The fields 'z', 'i' and 'p' correspond to the
%   nominal p.u. complex power (at 1 p.u. voltage magnitude) of the constant
%   impedance, constant current, and constant power portions, respectively of
%   the ZIP load model.
%
%   Example:
%       Sd = makeSdzip(baseMVA, bus, mpopt);

%   MATPOWER
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Shrirang Abhyankar
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

if nargin < 3
    mpopt = [];
end
if ~isempty(mpopt) && ~isempty(mpopt.exp.sys_wide_zip_loads.pw)
    if any(size(mpopt.exp.sys_wide_zip_loads.pw) ~= [1 3])
        error('makeSdzip: ''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
    end
    if abs(sum(mpopt.exp.sys_wide_zip_loads.pw) - 1) > eps
        error('makeSdzip: elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
    end
    pw = mpopt.exp.sys_wide_zip_loads.pw;
else
    pw = [1 0 0];
end
if ~isempty(mpopt) && ~isempty(mpopt.exp.sys_wide_zip_loads.qw)
    if any(size(mpopt.exp.sys_wide_zip_loads.qw) ~= [1 3])
        error('makeSdzip: ''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
    end
    if abs(sum(mpopt.exp.sys_wide_zip_loads.qw) - 1) > eps
        error('makeSdzip: elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
    end
    qw = mpopt.exp.sys_wide_zip_loads.qw;
else
    qw = pw;
end

Sd.z = (bus(:, PD) * pw(3)  + 1j * bus(:, QD) * qw(3)) / baseMVA;
Sd.i = (bus(:, PD) * pw(2)  + 1j * bus(:, QD) * qw(2)) / baseMVA;
Sd.p = (bus(:, PD) * pw(1)  + 1j * bus(:, QD) * qw(1)) / baseMVA;

function [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V)
%DSBUS_DV   Computes partial derivatives of power injection w.r.t. voltage.
%   [DSBUS_DVM, DSBUS_DVA] = DSBUS_DV(YBUS, V) returns two matrices containing
%   partial derivatives of the complex bus power injections w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If YBUS is a
%   sparse matrix, the return values will be also. The following explains
%   the expressions used to form the matrices:
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%
%   Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%       dI/dVm = Ybus * dV/dVm = Ybus * diag(V./abs(V))
%
%   Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%       dI/dVa = Ybus * dV/dVa = Ybus * j * diag(V)
%
%   Partials of S w.r.t. voltage magnitudes
%       dS/dVm = diag(V) * conj(dI/dVm) + diag(conj(Ibus)) * dV/dVm
%              = diag(V) * conj(Ybus * diag(V./abs(V)))
%                                       + conj(diag(Ibus)) * diag(V./abs(V))
%
%   Partials of S w.r.t. voltage angles
%       dS/dVa = diag(V) * conj(dI/dVa) + diag(conj(Ibus)) * dV/dVa
%              = diag(V) * conj(Ybus * j * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = -j * diag(V) * conj(Ybus * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = j * diag(V) * conj(diag(Ibus) - Ybus * diag(V))
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

n = length(V);
Ibus = Ybus * V;

if issparse(Ybus)           %% sparse version (if Ybus is sparse)
    diagV       = sparse(1:n, 1:n, V, n, n);
    diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
    diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
else                        %% dense version
    diagV       = diag(V);
    diagIbus    = diag(Ibus);
    diagVnorm   = diag(V./abs(V));
end

dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
dSbus_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV);


function [Pd, Qd] = total_load(bus, gen, load_zone, opt, mpopt)
%TOTAL_LOAD Returns vector of total load in each load zone.
%   PD = TOTAL_LOAD(MPC)
%   PD = TOTAL_LOAD(MPC, LOAD_ZONE)
%   PD = TOTAL_LOAD(MPC, LOAD_ZONE, OPT)
%   PD = TOTAL_LOAD(MPC, LOAD_ZONE, OPT, MPOPT)
%   PD = TOTAL_LOAD(BUS)
%   PD = TOTAL_LOAD(BUS, GEN)
%   PD = TOTAL_LOAD(BUS, GEN, LOAD_ZONE)
%   PD = TOTAL_LOAD(BUS, GEN, LOAD_ZONE, OPT)
%   PD = TOTAL_LOAD(BUS, GEN, LOAD_ZONE, OPT, MPOPT)
%   [PD, QD] = TOTAL_LOAD(...) returns both active and reative power
%   demand for each zone.
%
%   MPC - standard MATPOWER case struct
%
%   BUS - standard BUS matrix with nb rows, where the fixed active
%       and reactive loads are specified in columns PD and QD
%
%   GEN - (optional) standard GEN matrix with ng rows, where the
%       dispatchable loads are specified by columns PG, QG, PMIN,
%       QMIN and QMAX (in rows for which ISLOAD(GEN) returns true).
%       If GEN is empty, it assumes there are no dispatchable loads.
%
%   LOAD_ZONE - (optional) nb element vector where the value of
%       each element is either zero or the index of the load zone
%       to which the corresponding bus belongs. If LOAD_ZONE(b) = k
%       then the loads at bus b will added to the values of PD(k) and
%       QD(k). If LOAD_ZONE is empty, the default is defined as the areas
%       specified in the BUS matrix, i.e. LOAD_ZONE = BUS(:, BUS_AREA)
%       and load will have dimension = MAX(BUS(:, BUS_AREA)). LOAD_ZONE
%       can also take the following string values:
%           'all'  - use a single zone for the entire system (return scalar)
%           'area' - use LOAD_ZONE = BUS(:, BUS_AREA), same as default
%           'bus'  - use a different zone for each bus (i.e. to compute
%               final values of bus-wise loads, including voltage dependent
%               fixed loads and or dispatchable loads)
%
%   OPT - (optional) option struct, with the following fields:
%           'type'  -  string specifying types of loads to include, default
%                      is 'BOTH' if GEN is provided, otherwise 'FIXED'
%               'FIXED'        : sum only fixed loads
%               'DISPATCHABLE' : sum only dispatchable loads
%               'BOTH'         : sum both fixed and dispatchable loads
%           'nominal' -  1 : use nominal load for dispatchable loads
%                        0 : (default) use actual realized load for
%                             dispatchable loads
%
%       For backward compatibility with MATPOWER 4.x, OPT can also
%       take the form of a string, with the same options as OPT.type above.
%       In this case, again for backward compatibility, it is the "nominal"
%       load that is computed for dispatchable loads, not the actual
%       realized load. Using a string for OPT is deprecated and
%       will be removed in a future version.
%
%   MPOPT - (optional) MATPOWER options struct, which may specify
%       a voltage dependent (ZIP) load model for fixed loads
%
%   Examples:
%       Return the total active load for each area as defined in BUS_AREA.
%
%       Pd = total_load(bus);
%
%       Return total active and reactive load, fixed and dispatchable, for
%       entire system.
%
%       [Pd, Qd] = total_load(bus, gen, 'all');
%
%       Return the total of the nominal dispatchable loads at buses 10-20.
%
%       load_zone = zeros(nb, 1);
%       load_zone(10:20) = 1;
%       opt = struct('type', 'DISPATCHABLE', 'nominal', 1);
%       Pd = total_load(mpc, load_zone, opt)
%
%   See also SCALE_LOAD.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%% purposely being backward compatible with older MATPOWER
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%%-----  process inputs  -----
if nargin < 5
    mpopt = [];
    if nargin < 4
        opt = [];
        if nargin < 3
            load_zone = [];
            if nargin < 2
                gen = [];
            end
        end
    end
end
if isstruct(bus)    %% shift input arguments
    mpopt = opt;
    opt = load_zone;
    load_zone = gen;
    mpc = bus;
    bus = mpc.bus;
    gen = mpc.gen;
end

nb = size(bus, 1);          %% number of buses

%% default options
if ischar(opt)      %% convert old WHICH_TYPE string option to struct
    opt = struct('type', opt, 'nominal', 1);
else
    if ~isfield(opt, 'type')
        if isempty(gen)
            opt.type = 'FIXED';
        else
            opt.type = 'BOTH';
        end
    end
    if ~isfield(opt, 'nominal')
        opt.nominal = 0;
    end
end
switch upper(opt.type(1))
    case {'F', 'D', 'B'}
        %% OK
    otherwise
        error('total_load: OPT.type should be ''FIXED'', ''DISPATCHABLE'' or ''BOTH''');
end
want_Q      = (nargout > 1);
want_fixed  = (opt.type(1) == 'B' || opt.type(1) == 'F');
want_disp   = (opt.type(1) == 'B' || opt.type(1) == 'D');

%% initialize load_zone
if ischar(load_zone)
    if strcmp(lower(load_zone), 'bus')
        load_zone = (1:nb)';            %% each bus is its own zone
    elseif strcmp(lower(load_zone), 'all')
        load_zone = ones(nb, 1);        %% make a single zone of all buses
    elseif strcmp(lower(load_zone), 'area')
        load_zone = bus(:, BUS_AREA);   %% use areas defined in bus data as zones
    end
elseif isempty(load_zone)
    load_zone = bus(:, BUS_AREA);   %% use areas defined in bus data as zones
end
nz = max(load_zone);    %% number of load zones

%% fixed load at each bus, & initialize dispatchable
if want_fixed
    Sd = makeSdzip(1, bus, mpopt);
    Vm = bus(:, VM);
    Sbusd = Sd.p + Sd.i .* Vm + Sd.z .* Vm.^2;
    Pdf = real(Sbusd);      %% real power
    if want_Q
        Qdf = imag(Sbusd);  %% reactive power
    end
else
    Pdf = zeros(nb, 1);     %% real power
    if want_Q
        Qdf = zeros(nb, 1); %% reactive power
    end
end

%% dispatchable load at each bus 
if want_disp            %% need dispatchable
    ng = size(gen, 1);
    is_ld = isload(gen) & gen(:, GEN_STATUS) > 0;
    ld = find(is_ld);

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:nb)';

    Cld = sparse(e2i(gen(:, GEN_BUS)), (1:ng)', is_ld, nb, ng);
    if opt.nominal      %% use nominal power
        Pdd = -Cld * gen(:, PMIN);      %% real power
        if want_Q
            Q = zeros(ng, 1);
            Q(ld) = (gen(ld, QMIN) == 0) .* gen(ld, QMAX) + ...
                    (gen(ld, QMAX) == 0) .* gen(ld, QMIN);
            Qdd = -Cld * Q;             %% reactive power
        end
    else                %% use realized actual power dispatch
        Pdd = -Cld * gen(:, PG);        %% real power
        if want_Q
            Qdd = -Cld * gen(:, QG);    %% reactive power
        end
    end
else
    Pdd = zeros(nb, 1);
    if want_Q
        Qdd = zeros(nb, 1);
    end
end

%% compute load sums
if nz == nb && all(load_zone == (1:nb)');   %% individual buses
    Pd = (Pdf + Pdd) .* (bus(:, BUS_TYPE) ~= NONE);
    if want_Q
        Qd = (Qdf + Qdd) .* (bus(:, BUS_TYPE) ~= NONE);
    end
else
    Pd = zeros(nz, 1);
    if want_Q
        Qd = zeros(nz, 1);
    end
    for k = 1:nz
        idx = find( load_zone == k & bus(:, BUS_TYPE) ~= NONE);
        Pd(k) = sum(Pdf(idx)) + sum(Pdd(idx));
        if want_Q
            Qd(k) = sum(Qdf(idx)) + sum(Qdd(idx));
        end
    end
end

function [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas)
%INT2EXT   Converts internal to external bus numbering.
%
%   This function has two forms, (1) the old form that operates on
%   and returns individual matrices and (2) the new form that operates
%   on and returns an entire MATPOWER case struct.
%
%   1.  [BUS, GEN, BRANCH, AREAS] = INT2EXT(I2E, BUS, GEN, BRANCH, AREAS)
%       [BUS, GEN, BRANCH] = INT2EXT(I2E, BUS, GEN, BRANCH)
%
%   Converts from the consecutive internal bus numbers back to the originals
%   using the mapping provided by the I2E vector returned from EXT2INT,
%   where EXTERNAL_BUS_NUMBER = I2E(INTERNAL_BUS_NUMBER).
%   AREAS is completely ignored and is only included here for backward
%   compatibility of the API.
%
%   Examples:
%       [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
%       [bus, gen, branch] = int2ext(i2e, bus, gen, branch);
%
%   2.  MPC = INT2EXT(MPC)
%
%   If the input is a single MATPOWER case struct, then it restores all
%   buses, generators and branches that were removed because of being
%   isolated or off-line, and reverts to the original generator ordering
%   and original bus numbering. This requires that the 'order' field
%   created by EXT2INT be in place.
%
%   Example:
%       mpc = int2ext(mpc);
%
%   See also EXT2INT, I2E_FIELD, I2E_DATA.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isstruct(i2e)
    mpc = i2e;
    if nargin == 1
        if ~isfield(mpc, 'order')
            error('int2ext: mpc does not have the ''order'' field required for conversion back to external numbering.');
        end
        o = mpc.order;

        if o.state == 'i'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I] = idx_bus;
            GEN_BUS = idx_gen;
            [F_BUS, T_BUS] = idx_brch;

            %% execute userfcn callbacks for 'int2ext' stage
            if isfield(mpc, 'userfcn')
                mpc = run_userfcn(mpc.userfcn, 'int2ext', mpc);
            end
            
            %% convert back "extra" fields
            if isfield(mpc, 'gencost')
                ordering = {'gen'};         %% Pg cost only
                if size(mpc.gencost, 1) == 2*size(mpc.gen, 1)
                    ordering{2} = 'gen';    %% include Qg cost
                end
                mpc = i2e_field(mpc, 'gencost', ordering);
            end
            %% assume A and N are "read-only"
            %% (otherwise need to convert back, using i2e_field() which
            %% requires knowing if they are sized for AC or DC)
            if isfield(mpc, 'A')
                o.int.A = mpc.A;
                mpc.A = o.ext.A;
            end
            if isfield(mpc, 'N')
                o.int.N = mpc.N;
                mpc.N = o.ext.N;
            end

            %% save data matrices with internal ordering & restore originals
            o.int.bus    = mpc.bus;
            o.int.branch = mpc.branch;
            o.int.gen    = mpc.gen;
            mpc.bus     = o.ext.bus;
            mpc.branch  = o.ext.branch;
            mpc.gen     = o.ext.gen;

            %% zero pad data matrices on right if necessary
            nci = size(o.int.bus, 2);
            [nr, nc] = size(mpc.bus);
            if nc < nci
                mpc.bus = [mpc.bus zeros(nr, nci-nc)];
            end
            nci = size(o.int.branch, 2);
            [nr, nc] = size(mpc.branch);
            if nc < nci
                mpc.branch = [mpc.branch zeros(nr, nci-nc)];
            end
            nci = size(o.int.gen, 2);
            [nr, nc] = size(mpc.gen);
            if nc < nci
                mpc.gen = [mpc.gen zeros(nr, nci-nc)];
            end

            %% update data (in bus, branch, and gen only)
            mpc.bus(o.bus.status.on, :)       = o.int.bus;
            mpc.branch(o.branch.status.on, :) = o.int.branch;
            mpc.gen(o.gen.status.on, :)       = o.int.gen(o.gen.i2e, :);

            %% revert to original bus numbers
            mpc.bus(o.bus.status.on, BUS_I) = ...
                    o.bus.i2e( mpc.bus(o.bus.status.on, BUS_I) );
            mpc.branch(o.branch.status.on, F_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, F_BUS) );
            mpc.branch(o.branch.status.on, T_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, T_BUS) );
            mpc.gen(o.gen.status.on, GEN_BUS) = ...
                    o.bus.i2e( mpc.gen(o.gen.status.on, GEN_BUS) );

            if isfield(o, 'ext')
                o = rmfield(o, 'ext');
            end
            o.state = 'e';
            mpc.order = o;
        else
            error('int2ext: mpc claims it is already using external numbering.');
        end

        bus = mpc;
    else                    %% convert extra data
        if ischar(bus) || iscell(bus)   %% field
            warning('Calls of the form MPC = INT2EXT(MPC, ''FIELD_NAME'', ...) have been deprecated. Please replace INT2EXT with I2E_FIELD.');
            if nargin > 3
                dim = branch;
            else
                dim = 1;
            end
            bus = i2e_field(mpc, bus, gen, dim);
        else                            %% value
            warning('Calls of the form VAL = INT2EXT(MPC, VAL, ...) have been deprecated. Please replace INT2EXT with I2E_DATA.');
            if nargin > 4
                dim = areas;
            else
                dim = 1;
            end
            bus = i2e_data(mpc, bus, gen, branch, dim);
        end
    end
else            %% old form
    %% define names for columns to data matrices
    [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    [GEN_BUS] = idx_gen;
    [F_BUS, T_BUS] = idx_brch;

    bus(:, BUS_I)               = i2e( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = i2e( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = i2e( branch(:, F_BUS)         );
    branch(:, T_BUS)            = i2e( branch(:, T_BUS)         );
end

function newval = i2e_data(mpc, val, oldval, ordering, dim)
%I2E_DATA   Converts data from internal to external bus numbering.
%
%   VAL = I2E_DATA(MPC, VAL, OLDVAL, ORDERING)
%   VAL = I2E_DATA(MPC, VAL, OLDVAL, ORDERING, DIM)
%
%   For a case struct using internal indexing, this function can be
%   used to convert other data structures as well by passing in 3 or 4
%   extra parameters in addition to the case struct. If the value passed
%   in the 2nd argument (VAL) is a column vector or cell array, it will
%   be converted according to the ordering specified by the 4th argument
%   (ORDERING, described below). If VAL is an n-dimensional matrix or
%   cell array, then the optional 5th argument (DIM, default = 1) can be
%   used to specify which dimension to reorder. The 3rd argument (OLDVAL)
%   is used to initialize the return value before converting VAL to
%   external indexing. In particular, any data corresponding to off-line
%   gens or branches or isolated buses or any connected gens or branches
%   will be taken from OLDVAL, with VAL supplying the rest of the
%   returned data.
%
%   The ORDERING argument is used to indicate whether the data
%   corresponds to bus-, gen- or branch-ordered data. It can be one
%   of the following three strings: 'bus', 'gen' or 'branch'. For
%   data structures with multiple blocks of data, ordered by bus,
%   gen or branch, they can be converted with a single call by
%   specifying ORDERING as a cell array of strings.
%
%   Any extra elements, rows, columns, etc. beyond those indicated
%   in ORDERING, are not disturbed.
%
%   Examples:
%       A_ext = i2e_data(mpc, A_int, A_orig, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       internal to external ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_ext = i2e_data(mpc, gencost_int, gencost_orig, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%   See also E2I_DATA, I2E_FIELD, INT2EXT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if ~isfield(mpc, 'order')
    error('i2e_data: mpc does not have the ''order'' field required for conversion back to external numbering.');
end
o = mpc.order;
if o.state ~= 'i'
    error('i2e_data: mpc does not appear to be in internal order');
end
if nargin < 5
    dim = 1;
end
if ischar(ordering)         %% single set
    if strcmp(ordering, 'gen')
        v = get_reorder(val, o.(ordering).i2e, dim);
    else
        v = val;
    end
    newval = set_reorder(oldval, v, o.(ordering).status.on, dim);
else                        %% multiple sets
    be = 0;  %% base, external indexing
    bi = 0;  %% base, internal indexing
    for k = 1:length(ordering)
        ne = size(o.ext.(ordering{k}), 1);
        ni = size(mpc.(ordering{k}), 1);
        v = get_reorder(val, bi+(1:ni), dim);
        oldv = get_reorder(oldval, be+(1:ne), dim);
        new_v{k} = i2e_data(mpc, v, oldv, ordering{k}, dim);
        be = be + ne;
        bi = bi + ni;
    end
    ni = size(val, dim);
    if ni > bi              %% the rest
        v = get_reorder(val, bi+1:ni, dim);
        new_v{length(new_v)+1} = v;
    end
    newval = cat(dim, new_v{:});
end

function mpc = i2e_field(mpc, field, ordering, dim)
%I2E_FIELD   Converts fields of MPC from internal to external bus numbering.
%
%   MPC = I2E_FIELD(MPC, FIELD, ORDERING)
%   MPC = I2E_FIELD(MPC, FIELD, ORDERING, DIM)
%
%   For a case struct using internal indexing, this function can be
%   used to convert other data structures as well by passing in 2 or 3
%   extra parameters in addition to the case struct.
%
%   The 2nd argument is a string or cell array of strings, specifying
%   a field in the case struct whose value should be converted by
%   a corresponding call to I2E_DATA. The field can contain either a
%   numeric or a cell array. The corresponding OLDVAL is taken from
%   where it was stored by EXT2INT in MPC.ORDER.EXT and the updated
%   case struct is returned. If FIELD is a cell array of strings,
%   they specify nested fields.
%
%   The 3rd and optional 4th arguments are simply passed along to
%   the call to I2E_DATA.
%
%   Examples:
%       mpc = i2e_field(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match external generator
%       ordering.
%
%       mpc = i2e_field(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match external
%       generator ordering.
%
%   See also E2I_FIELD, I2E_DATA, INT2EXT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    dim = 1;
end
if ischar(field)
    mpc.order.int.(field) = mpc.(field);
    mpc.(field) = i2e_data(mpc, mpc.(field), ...
                    mpc.order.ext.(field), ordering, dim);
else    %% iscell(field)
    for k = 1:length(field)
        s(k).type = '.';
        s(k).subs = field{k};
    end
    if ~isfield(mpc.order, 'int')
        mpc.order.int = [];
    end
    mpc.order.int = subsasgn(mpc.order.int, s, subsref(mpc, s));
    mpc = subsasgn(mpc, s, i2e_data(mpc, subsref(mpc, s), ...
        subsref(mpc.order.ext, s), ordering, dim));
end

function A = set_reorder(A, B, idx, dim)
%SET_REORDER Assigns B to A with one of the dimensions of A indexed.
%
%   A = SET_REORDER(A, B, IDX, DIM)
%
%   Returns A after doing A(:, ..., :, IDX, :, ..., :) = B
%   where DIM determines in which dimension to place the IDX.
%
%   If any dimension of B is smaller than the corresponding dimension
%   of A, the "extra" elements in A are untouched. If any dimension of
%   B is larger than the corresponding dimension of A, then A is padded
%   with zeros (if numeric) or empty matrices (if cell array) before
%   performing the assignment.
%
%   See also GET_REORDER.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check dimensions
ndim = ndims(A);
sA = size(A);
sB = size(B);
d = (1:length(sA));
d(dim) = [];        %% indices of all dimensions other than DIM

%% pad A with zeros (numeric) or empty matrices (cell), if necessary
s.subs = cell(1, ndim);
if any(sA(d) < sB(d))
    s.subs = num2cell(max(sA, sB));
    if iscell(A)
        s.type = '{}';
        A = subsasgn(A, s, []);
    else
        s.type = '()';
        A = subsasgn(A, s, 0);
    end
end

%% set up indexing
s.type = '()';
for k = 1:ndim
    if k == dim
        s.subs{k} = idx;
    else
        if sA(k) == sB(k)
            s.subs{k} = ':';        %% indexes of all elements in this dimension
        else    %% sA(k) > sB(k)
            s.subs{k} = (1:sB(k));  %% limit indexes to smaller size of B
        end
    end
end

%% do the assignment
A = subsasgn(A, s, B);


function printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt)
%PRINTPF   Prints power flow results.
%   PRINTPF(RESULTS, FD, MPOPT)
%   PRINTPF(BASEMVA, BUS, GEN, BRANCH, F, SUCCESS, ET, FD, MPOPT)
%
%   Prints power flow and optimal power flow results to FD (a file
%   descriptor which defaults to STDOUT), with the details of what
%   gets printed controlled by the optional MPOPT argument, which is a
%   MATPOWER options struct (see MPOPTION for details).
%
%   The data can either be supplied in a single RESULTS struct, or
%   in the individual arguments: BASEMVA, BUS, GEN, BRANCH, F, SUCCESS
%   and ET, where F is the OPF objective function value, SUCCESS is
%   true if the solution converged and false otherwise, and ET is the
%   elapsed time for the computation in seconds. If F is given, it is
%   assumed that the output is from an OPF run, otherwise it is assumed
%   to be a simple power flow run.
%
%   Examples:
%       mpopt = mpoptions('out.gen', 1, 'out.bus', 0, 'out.branch', 0);
%       [fd, msg] = fopen(fname, 'at');
%       results = runopf(mpc);
%       printpf(results);
%       printpf(results, fd);
%       printpf(results, fd, mpopt);
%       printpf(baseMVA, bus, gen, branch, f, success, et);
%       printpf(baseMVA, bus, gen, branch, f, success, et, fd);
%       printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt);
%       fclose(fd);

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default arguments
if isstruct(baseMVA)
    have_results_struct = 1;
    results = baseMVA;
    if nargin < 3 || isempty(gen)
        mpopt = mpoption;   %% use default options
    else
        mpopt = gen;
    end
    if mpopt.out.all == 0
        return;         %% nothin' to see here, bail out now
    end
    if nargin < 2 || isempty(bus)
        fd = 1;         %% print to stdio by default
    else
        fd = bus;
    end
    [baseMVA, bus, gen, branch, success, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, ...
            results.success, results.et);
    if isfield(results, 'f') && ~isempty(results.f)
        f = results.f;
    else
        f = [];
    end
else
    have_results_struct = 0;
    if nargin < 9
        mpopt = mpoption;   %% use default options
        if nargin < 8
            fd = 1;         %% print to stdio by default
        end
    end
    if mpopt.out.all == 0
        return;         %% nothin' to see here, bail out now
    end
end
isOPF = ~isempty(f);    %% FALSE -> only simple PF data, TRUE -> OPF data

%% options
isDC            = strcmp(upper(mpopt.model), 'DC');

SUPPRESS        = mpopt.out.suppress_detail;
if SUPPRESS == -1
    if size(bus, 1) > 500
        SUPPRESS = 1;
    else
        SUPPRESS = 0;
    end
end
OUT_ALL         = mpopt.out.all;
OUT_FORCE       = mpopt.out.force;
OUT_ANY         = OUT_ALL == 1;     %% set to true if any pretty output is to be generated
OUT_SYS_SUM     = OUT_ALL == 1 || (OUT_ALL == -1 && mpopt.out.sys_sum);
OUT_AREA_SUM    = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.area_sum);
OUT_BUS         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.bus);
OUT_BRANCH      = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.branch);
OUT_GEN         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && mpopt.out.gen);
OUT_ANY         = OUT_ANY || (OUT_ALL == -1 && ...
                    (OUT_SYS_SUM || OUT_AREA_SUM || OUT_BUS || ...
                    OUT_BRANCH || OUT_GEN));
if OUT_ALL == -1
    OUT_ALL_LIM = ~SUPPRESS * mpopt.out.lim.all;
elseif OUT_ALL == 1
    OUT_ALL_LIM = 2;
else
    OUT_ALL_LIM = 0;
end
OUT_ANY         = OUT_ANY || OUT_ALL_LIM >= 1;
if OUT_ALL_LIM == -1
    OUT_V_LIM       = ~SUPPRESS * mpopt.out.lim.v;
    OUT_LINE_LIM    = ~SUPPRESS * mpopt.out.lim.line;
    OUT_PG_LIM      = ~SUPPRESS * mpopt.out.lim.pg;
    OUT_QG_LIM      = ~SUPPRESS * mpopt.out.lim.qg;
else
    OUT_V_LIM       = OUT_ALL_LIM;
    OUT_LINE_LIM    = OUT_ALL_LIM;
    OUT_PG_LIM      = OUT_ALL_LIM;
    OUT_QG_LIM      = OUT_ALL_LIM;
end
OUT_ANY         = OUT_ANY || (OUT_ALL_LIM == -1 && (OUT_V_LIM || OUT_LINE_LIM || OUT_PG_LIM || OUT_QG_LIM));

%%----- print the stuff -----
if OUT_ANY
    ptol = 1e-4;        %% tolerance for displaying shadow prices
    if isOPF && ~isDC && strcmp(upper(mpopt.opf.ac.solver), 'SDPOPF')
        isSDP = 1;
        ptol = 0.1;     %% tolerance for displaying shadow prices
        if have_results_struct && isfield(results, 'mineigratio') && ~isempty(results.mineigratio)
            mineigratio = results.mineigratio;
        else
            mineigratio = [];
        end
        if have_results_struct && isfield(results, 'zero_eval') && ~isempty(results.zero_eval)
            zero_eval = results.zero_eval;
        else
            zero_eval = [];
        end
    else
        isSDP = 0;
    end

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';

    %% sizes of things
    nb = size(bus, 1);      %% number of buses
    nl = size(branch, 1);   %% number of branches
    ng = size(gen, 1);      %% number of generators

    %% zero out some data to make printout consistent for DC case
    if isDC
        bus(:, [QD, BS])            = zeros(nb, 2);
        gen(:, [QG, QMAX, QMIN])    = zeros(ng, 3);
        branch(:, [BR_R, BR_B])     = zeros(nl, 2);
    end

    %% parameters
    ties = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= bus(e2i(branch(:, T_BUS)), BUS_AREA));
                            %% area inter-ties
    xfmr = find(branch(:, TAP));                    %% indices of transformers
    nzld = find((bus(:, PD) | bus(:, QD)) & bus(:, BUS_TYPE) ~= NONE);
    sorted_areas = sort(bus(:, BUS_AREA));
    s_areas = sorted_areas([1; find(diff(sorted_areas))+1]);    %% area numbers
    nzsh = find((bus(:, GS) | bus(:, BS)) & bus(:, BUS_TYPE) ~= NONE);
    allg = find( ~isload(gen) );
    alld = find(  isload(gen) );
    ong  = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );
    onld = find( gen(:, GEN_STATUS) > 0 &  isload(gen) );
    V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
    opt = struct('type', 'FIXED');
    [Pdf, Qdf] = total_load(bus, gen, 'bus', struct('type', 'FIXED'), mpopt);
    [Pdd, Qdd] = total_load(bus, gen, 'bus', struct('type', 'DISPATCHABLE'), mpopt);
    if isDC
        loss = zeros(nl, 1);
        fchg = loss;
        tchg = loss;
    else
        [loss, fchg, tchg] = get_losses(baseMVA, bus, branch);
    end

    %% convergence & elapsed time
    if success
        if isSDP
            fprintf(fd, '\nSolution satisfies rank and consistency conditions, %.2f seconds.\nmineigratio = %0.5g, zero_eval = %0.5g', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\nConverged in %.2f seconds', et);
        end
    else
        if isSDP
            fprintf(fd, '\n>>>>>  Solution does NOT satisfy rank and/or consistency conditions (%.2f seconds).  <<<<<\nmineigratio = %0.5g, zero_eval = %0.5g\n', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\n>>>>>  Did NOT converge (%.2f seconds)  <<<<<\n', et);
        end
    end

    %% objective function value
    if isOPF && (success || OUT_FORCE)
        fprintf(fd, '\nObjective Function Value = %.2f $/hr', f);
    end
end
if OUT_SYS_SUM && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     System Summary                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n\nHow many?                How much?              P (MW)            Q (MVAr)');
    fprintf(fd, '\n---------------------    -------------------  -------------  -----------------');
    fprintf(fd, '\nBuses         %6d     Total Gen Capacity   %7.1f       %7.1f to %.1f', nb, sum(gen(allg, PMAX)), sum(gen(allg, QMIN)), sum(gen(allg, QMAX)));
    fprintf(fd, '\nGenerators     %5d     On-line Capacity     %7.1f       %7.1f to %.1f', length(allg), sum(gen(ong, PMAX)), sum(gen(ong, QMIN)), sum(gen(ong, QMAX)));
    fprintf(fd, '\nCommitted Gens %5d     Generation (actual)  %7.1f           %7.1f', length(ong), sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\nLoads          %5d     Load                 %7.1f           %7.1f', length(nzld)+length(onld), sum(Pdf(nzld))-sum(gen(onld, PG)), sum(Qdf(nzld))-sum(gen(onld, QG)));
    fprintf(fd, '\n  Fixed        %5d       Fixed              %7.1f           %7.1f', length(nzld), sum(Pdf(nzld)), sum(Qdf(nzld)));
    fprintf(fd, '\n  Dispatchable %5d       Dispatchable       %7.1f of %-7.1f%7.1f', length(onld), -sum(gen(onld, PG)), -sum(gen(onld, PMIN)), -sum(gen(onld, QG)));
    fprintf(fd, '\nShunts         %5d     Shunt (inj)          %7.1f           %7.1f', length(nzsh), ...
        -sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, GS)), sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, BS)) );
    fprintf(fd, '\nBranches       %5d     Losses (I^2 * Z)     %8.2f          %8.2f', nl, sum(real(loss)), sum(imag(loss)) );
    fprintf(fd, '\nTransformers   %5d     Branch Charging (inj)     -            %7.1f', length(xfmr), sum(fchg) + sum(tchg) );
    fprintf(fd, '\nInter-ties     %5d     Total Inter-tie Flow %7.1f           %7.1f', length(ties), sum(abs(branch(ties, PF)-branch(ties, PT))) / 2, sum(abs(branch(ties, QF)-branch(ties, QT))) / 2);
    fprintf(fd, '\nAreas          %5d', length(s_areas));
    fprintf(fd, '\n');
    fprintf(fd, '\n                          Minimum                      Maximum');
    fprintf(fd, '\n                 -------------------------  --------------------------------');
    [minv, mini] = min(bus(:, VM));
    [maxv, maxi] = max(bus(:, VM));
    fprintf(fd, '\nVoltage Magnitude %7.3f p.u. @ bus %-4d     %7.3f p.u. @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    [minv, mini] = min(bus(:, VA));
    [maxv, maxi] = max(bus(:, VA));
    fprintf(fd, '\nVoltage Angle   %8.2f deg   @ bus %-4d   %8.2f deg   @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    if ~isDC
        [maxv, maxi] = max(real(loss));
        fprintf(fd, '\nP Losses (I^2*R)             -              %8.2f MW    @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
        [maxv, maxi] = max(imag(loss));
        fprintf(fd, '\nQ Losses (I^2*X)             -              %8.2f MVAr  @ line %d-%d', maxv, branch(maxi, F_BUS), branch(maxi, T_BUS));
    end
    if isOPF
        [minv, mini] = min(bus(:, LAM_P));
        [maxv, maxi] = max(bus(:, LAM_P));
        fprintf(fd, '\nLambda P        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
        [minv, mini] = min(bus(:, LAM_Q));
        [maxv, maxi] = max(bus(:, LAM_Q));
        fprintf(fd, '\nLambda Q        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d', minv, bus(mini, BUS_I), maxv, bus(maxi, BUS_I));
    end
    fprintf(fd, '\n');
end

if OUT_AREA_SUM && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Area Summary                                                             |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nArea  # of      # of Gens        # of Loads         # of    # of   # of   # of');
    fprintf(fd, '\n Num  Buses   Total  Online   Total  Fixed  Disp    Shunt   Brchs  Xfmrs   Ties');
    fprintf(fd, '\n----  -----   -----  ------   -----  -----  -----   -----   -----  -----  -----');
    for i=1:length(s_areas)
        a = s_areas(i);
        ib = find(bus(:, BUS_AREA) == a);
        ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & ~isload(gen));
        igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & ~isload(gen));
        ildon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & isload(gen));
        inzld = find(bus(:, BUS_AREA) == a & (Pdf | Qdf) & bus(:, BUS_TYPE) ~= NONE);
        inzsh = find(bus(:, BUS_AREA) == a & (bus(:, GS) | bus(:, BS)) & bus(:, BUS_TYPE) ~= NONE);
        ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a);
        out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a);
        if isempty(xfmr)
            nxfmr = 0;
        else
            nxfmr = length(find(bus(e2i(branch(xfmr, F_BUS)), BUS_AREA) == a & bus(e2i(branch(xfmr, T_BUS)), BUS_AREA) == a));
        end
        fprintf(fd, '\n%3d  %6d   %5d  %5d   %5d  %5d  %5d   %5d   %5d  %5d  %5d', ...
            a, length(ib), length(ig), length(igon), ...
            length(inzld)+length(ildon), length(inzld), length(ildon), ...
            length(inzsh), length(ibrch), nxfmr, length(in_tie)+length(out_tie));
    end
    fprintf(fd, '\n----  -----   -----  ------   -----  -----  -----   -----   -----  -----  -----');
    fprintf(fd, '\nTot: %6d   %5d  %5d   %5d  %5d  %5d   %5d   %5d  %5d  %5d', ...
        nb, length(allg), length(ong), length(nzld)+length(onld), ...
        length(nzld), length(onld), length(nzsh), nl, length(xfmr), length(ties));
    fprintf(fd, '\n');
    fprintf(fd, '\nArea      Total Gen Capacity           On-line Gen Capacity         Generation');
    fprintf(fd, '\n Num     MW           MVAr            MW           MVAr             MW    MVAr');
    fprintf(fd, '\n----   ------  ------------------   ------  ------------------    ------  ------');
    for i=1:length(s_areas)
        a = s_areas(i);
        ig = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & ~isload(gen));
        igon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & ~isload(gen));
        fprintf(fd, '\n%3d   %7.1f  %7.1f to %-7.1f  %7.1f  %7.1f to %-7.1f   %7.1f %7.1f', ...
            a, sum(gen(ig, PMAX)), sum(gen(ig, QMIN)), sum(gen(ig, QMAX)), ...
            sum(gen(igon, PMAX)), sum(gen(igon, QMIN)), sum(gen(igon, QMAX)), ...
            sum(gen(igon, PG)), sum(gen(igon, QG)) );
    end
    fprintf(fd, '\n----   ------  ------------------   ------  ------------------    ------  ------');
    fprintf(fd, '\nTot:  %7.1f  %7.1f to %-7.1f  %7.1f  %7.1f to %-7.1f   %7.1f %7.1f', ...
            sum(gen(allg, PMAX)), sum(gen(allg, QMIN)), sum(gen(allg, QMAX)), ...
            sum(gen(ong, PMAX)), sum(gen(ong, QMIN)), sum(gen(ong, QMAX)), ...
            sum(gen(ong, PG)), sum(gen(ong, QG)) );
    fprintf(fd, '\n');
    fprintf(fd, '\nArea    Disp Load Cap       Disp Load         Fixed Load        Total Load');
    fprintf(fd, '\n Num      MW     MVAr       MW     MVAr       MW     MVAr       MW     MVAr');
    fprintf(fd, '\n----    ------  ------    ------  ------    ------  ------    ------  ------');
    Qlim = (gen(:, QMIN) == 0) .* gen(:, QMAX) + (gen(:, QMAX) == 0) .* gen(:, QMIN);
    for i=1:length(s_areas)
        a = s_areas(i);
        ildon = find(bus(e2i(gen(:, GEN_BUS)), BUS_AREA) == a & gen(:, GEN_STATUS) > 0 & isload(gen));
        inzld = find(bus(:, BUS_AREA) == a & (Pdf | Qdf));
        fprintf(fd, '\n%3d    %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            a, -sum(gen(ildon, PMIN)), ...
            -sum(Qlim(ildon)), ...
            -sum(gen(ildon, PG)), -sum(gen(ildon, QG)), ...
            sum(Pdf(inzld)), sum(Qdf(inzld)), ...
            -sum(gen(ildon, PG)) + sum(Pdf(inzld)), ...
            -sum(gen(ildon, QG)) + sum(Qdf(inzld)) );
    end
    fprintf(fd, '\n----    ------  ------    ------  ------    ------  ------    ------  ------');
    fprintf(fd, '\nTot:   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f   %7.1f %7.1f', ...
            -sum(gen(onld, PMIN)), ...
            -sum(Qlim(onld)), ...
            -sum(gen(onld, PG)), -sum(gen(onld, QG)), ...
            sum(Pdf(nzld)), sum(Qdf(nzld)), ...
            -sum(gen(onld, PG)) + sum(Pdf(nzld)), ...
            -sum(gen(onld, QG)) + sum(Qdf(nzld)) );
    fprintf(fd, '\n');
    fprintf(fd, '\nArea      Shunt Inj        Branch      Series Losses      Net Export');
    fprintf(fd, '\n Num      MW     MVAr     Charging      MW     MVAr       MW     MVAr');
    fprintf(fd, '\n----    ------  ------    --------    ------  ------    ------  ------');
    for i=1:length(s_areas)
        a = s_areas(i);
        inzsh = find(bus(:, BUS_AREA) == a & (bus(:, GS) | bus(:, BS)));
        ibrch = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
        in_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) ~= a & bus(e2i(branch(:, T_BUS)), BUS_AREA) == a & branch(:, BR_STATUS));
        out_tie = find(bus(e2i(branch(:, F_BUS)), BUS_AREA) == a & bus(e2i(branch(:, T_BUS)), BUS_AREA) ~= a & branch(:, BR_STATUS));
        fprintf(fd, '\n%3d    %7.1f %7.1f    %7.1f    %7.2f %7.2f   %7.1f %7.1f', ...
            a, -sum(bus(inzsh, VM) .^ 2 .* bus(inzsh, GS)), ...
            sum(bus(inzsh, VM) .^ 2 .* bus(inzsh, BS)), ...
            sum(fchg(ibrch)) + sum(tchg(ibrch)) + sum(fchg(out_tie)) + sum(tchg(in_tie)), ...
            sum(real(loss(ibrch))) + sum(real(loss([in_tie; out_tie]))) / 2, ...
            sum(imag(loss(ibrch))) + sum(imag(loss([in_tie; out_tie]))) / 2, ...
            sum(branch(in_tie, PT))+sum(branch(out_tie, PF)) - sum(real(loss([in_tie; out_tie]))) / 2, ...
            sum(branch(in_tie, QT))+sum(branch(out_tie, QF)) - sum(imag(loss([in_tie; out_tie]))) / 2  );
    end
    fprintf(fd, '\n----    ------  ------    --------    ------  ------    ------  ------');
    fprintf(fd, '\nTot:   %7.1f %7.1f    %7.1f    %7.2f %7.2f       -       -', ...
        -sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, GS)), ...
        sum(bus(nzsh, VM) .^ 2 .* bus(nzsh, BS)), ...
        sum(fchg) + sum(tchg), sum(real(loss)), sum(imag(loss)) );
    fprintf(fd, '\n');
end

%% generator data
if OUT_GEN && (success || OUT_FORCE)
    if isOPF
        genlamP = bus(e2i(gen(:, GEN_BUS)), LAM_P);
        genlamQ = bus(e2i(gen(:, GEN_BUS)), LAM_Q);
    end
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Generator Data                                                           |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Gen   Bus   Status     Pg        Qg   ');
    if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
    fprintf(fd, '\n  #     #              (MW)     (MVAr) ');
    if isOPF, fprintf(fd, '     P         Q    '); end
    fprintf(fd, '\n----  -----  ------  --------  --------');
    if isOPF, fprintf(fd, '  --------  --------'); end
    for k = 1:length(allg)
        i = allg(k);
        fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
        if gen(i, GEN_STATUS) > 0 && (gen(i, PG) || gen(i, QG))
            fprintf(fd, '%10.2f%10.2f', gen(i, PG), gen(i, QG));
        else
            fprintf(fd, '       -         -  ');
        end
        if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
    end
    fprintf(fd, '\n                     --------  --------');
    fprintf(fd, '\n            Total: %9.2f%10.2f', sum(gen(ong, PG)), sum(gen(ong, QG)));
    fprintf(fd, '\n');
    if ~isempty(alld)
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Data                                                   |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n Gen   Bus   Status     Pd        Qd   ');
        if isOPF, fprintf(fd, '   Lambda ($/MVA-hr)'); end
        fprintf(fd, '\n  #     #              (MW)     (MVAr) ');
        if isOPF, fprintf(fd, '     P         Q    '); end
        fprintf(fd, '\n----  -----  ------  --------  --------');
        if isOPF, fprintf(fd, '  --------  --------'); end
        for k = 1:length(alld)
            i = alld(k);
            fprintf(fd, '\n%3d %6d     %2d ', i, gen(i, GEN_BUS), gen(i, GEN_STATUS));
            if gen(i, GEN_STATUS) > 0 && (gen(i, PG) || gen(i, QG))
                fprintf(fd, '%10.2f%10.2f', -gen(i, PG), -gen(i, QG));
            else
                fprintf(fd, '       -         -  ');
            end
            if isOPF, fprintf(fd, '%10.2f%10.2f', genlamP(i), genlamQ(i)); end
        end
        fprintf(fd, '\n                     --------  --------');
        fprintf(fd, '\n            Total: %9.2f%10.2f', -sum(gen(onld, PG)), -sum(gen(onld, QG)));
        fprintf(fd, '\n');
    end
end

%% bus data
if OUT_BUS && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Bus Data                                                                 |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Bus      Voltage          Generation             Load        ');
    if isOPF, fprintf(fd, '  Lambda($/MVA-hr)'); end
    fprintf(fd, '\n  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
    if isOPF, fprintf(fd, '     P        Q   '); end
    fprintf(fd, '\n----- ------- --------  --------  --------  --------  --------');
    if isOPF, fprintf(fd, '  -------  -------'); end
    for i = 1:nb
        fprintf(fd, '\n%5d%7.3f%9.3f', bus(i, [BUS_I, VM, VA]));
        if bus(i, BUS_TYPE) == REF
            fprintf(fd, '*');
        elseif bus(i, BUS_TYPE) == NONE
            fprintf(fd, 'x');
        else
            fprintf(fd, ' ');
        end
        g  = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    ~isload(gen));
        ld = find(gen(:, GEN_STATUS) > 0 & gen(:, GEN_BUS) == bus(i, BUS_I) & ...
                    isload(gen));
        if ~isempty(g)
            fprintf(fd, '%9.2f%10.2f', sum(gen(g, PG)), sum(gen(g, QG)));
        else
            fprintf(fd, '      -         -  ');
        end
        if Pdf(i) || Qdf(i) || ~isempty(ld)
            if ~isempty(ld)
                fprintf(fd, '%10.2f*%9.2f*', Pdf(i) - sum(gen(ld, PG)), ...
                                             Qdf(i) - sum(gen(ld, QG)));
            else
                fprintf(fd, '%10.2f%10.2f ', [ Pdf(i) Qdf(i) ]);
            end
        else
            fprintf(fd, '       -         -   ');
        end
        if isOPF
            fprintf(fd, '%9.3f', bus(i, LAM_P));
            if abs(bus(i, LAM_Q)) > ptol
                fprintf(fd, '%8.3f', bus(i, LAM_Q));
            else
                fprintf(fd, '     -');
            end
        end
    end
    fprintf(fd, '\n                        --------  --------  --------  --------');
    fprintf(fd, '\n               Total: %9.2f %9.2f %9.2f %9.2f', ...
        sum(gen(ong, PG)), sum(gen(ong, QG)), ...
        sum(Pdf(nzld)) - sum(gen(onld, PG)), ...
        sum(Qdf(nzld)) - sum(gen(onld, QG)));
    fprintf(fd, '\n');
end

%% branch data
if OUT_BRANCH && (success || OUT_FORCE)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Branch Data                                                              |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\nBrnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  ');
    fprintf(fd, '\n  #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)');
    fprintf(fd, '\n-----  -----  -----  --------  --------  --------  --------  --------  --------');
    fprintf(fd, '\n%4d%7d%7d%10.2f%10.2f%10.2f%10.2f%10.3f%10.2f', ...
            [   (1:nl)', branch(:, [F_BUS, T_BUS]), ...
                branch(:, [PF, QF]), branch(:, [PT, QT]), ...
                real(loss), imag(loss) ...
            ]');
    fprintf(fd, '\n                                                             --------  --------');
    fprintf(fd, '\n                                                    Total:%10.3f%10.2f', ...
            sum(real(loss)), sum(imag(loss)));
    fprintf(fd, '\n');
end

%%-----  constraint data  -----
if isOPF && (success || OUT_FORCE)
    ctol = mpopt.opf.violation; %% constraint violation tolerance
    %% voltage constraints
    if ~isDC && (OUT_V_LIM == 2 || (OUT_V_LIM == 1 && ...
                         (any(bus(:, VM) < bus(:, VMIN) + ctol) || ...
                          any(bus(:, VM) > bus(:, VMAX) - ctol) || ...
                          any(bus(:, MU_VMIN) > ptol) || ...
                          any(bus(:, MU_VMAX) > ptol))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Voltage Constraints                                                      |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBus #  Vmin mu    Vmin    |V|   Vmax    Vmax mu');
        fprintf(fd, '\n-----  --------   -----  -----  -----   --------');
        for i = 1:nb
            if OUT_V_LIM == 2 || (OUT_V_LIM == 1 && ...
                         (bus(i, VM) < bus(i, VMIN) + ctol || ...
                          bus(i, VM) > bus(i, VMAX) - ctol || ...
                          bus(i, MU_VMIN) > ptol || bus(i, MU_VMAX) > ptol))
                fprintf(fd, '\n%5d', bus(i, BUS_I));
                if bus(i, VM) < bus(i, VMIN) + ctol || bus(i, MU_VMIN) > ptol
                    fprintf(fd, '%10.3f', bus(i, MU_VMIN));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%8.3f%7.3f%7.3f', bus(i, [VMIN, VM, VMAX]));
                if bus(i, VM) > bus(i, VMAX) - ctol || bus(i, MU_VMAX) > ptol
                    fprintf(fd, '%10.3f', bus(i, MU_VMAX));
                else
                    fprintf(fd, '      -    ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% generator P constraints
    if OUT_PG_LIM == 2 || ...
            (OUT_PG_LIM == 1 && (any(gen(ong, PG) < gen(ong, PMIN) + ctol) || ...
                                any(gen(ong, PG) > gen(ong, PMAX) - ctol) || ...
                                any(gen(ong, MU_PMIN) > ptol) || ...
                                any(gen(ong, MU_PMAX) > ptol))) || ...
            (~isDC && (OUT_QG_LIM == 2 || ...
            (OUT_QG_LIM == 1 && (any(gen(ong, QG) < gen(ong, QMIN) + ctol) || ...
                                any(gen(ong, QG) > gen(ong, QMAX) - ctol) || ...
                                any(gen(ong, MU_QMIN) > ptol) || ...
                                any(gen(ong, MU_QMAX) > ptol)))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Generation Constraints                                                   |');
        fprintf(fd, '\n================================================================================');
    end
    if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                             (any(gen(ong, PG) < gen(ong, PMIN) + ctol) || ...
                              any(gen(ong, PG) > gen(ong, PMAX) - ctol) || ...
                              any(gen(ong, MU_PMIN) > ptol) || ...
                              any(gen(ong, MU_PMAX) > ptol)))
        fprintf(fd, '\n Gen   Bus                  Active Power Limits');
        fprintf(fd, '\n  #     #     Pmin mu     Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                        (gen(i, PG) < gen(i, PMIN) + ctol || ...
                         gen(i, PG) > gen(i, PMAX) - ctol || ...
                         gen(i, MU_PMIN) > ptol || gen(i, MU_PMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol || gen(i, MU_PMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol || gen(i, MU_PMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% generator Q constraints
    if ~isDC && (OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                             (any(gen(ong, QG) < gen(ong, QMIN) + ctol) || ...
                              any(gen(ong, QG) > gen(ong, QMAX) - ctol) || ...
                              any(gen(ong, MU_QMIN) > ptol) || ...
                              any(gen(ong, MU_QMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                 Reactive Power Limits');
        fprintf(fd, '\n  #     #     Qmin mu     Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(ong)
            i = ong(k);
            if OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                        (gen(i, QG) < gen(i, QMIN) + ctol || ...
                         gen(i, QG) > gen(i, QMAX) - ctol || ...
                         gen(i, MU_QMIN) > ptol || gen(i, MU_QMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol || gen(i, MU_QMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol || gen(i, MU_QMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% dispatchable load P constraints
    if ~isempty(onld) && (OUT_PG_LIM == 2 || ...
            (OUT_PG_LIM == 1 && (any(gen(onld, PG) < gen(onld, PMIN) + ctol) || ...
                                any(gen(onld, PG) > gen(onld, PMAX) - ctol) || ...
                                any(gen(onld, MU_PMIN) > ptol) || ...
                                any(gen(onld, MU_PMAX) > ptol))) || ...
            (~isDC && (OUT_QG_LIM == 2 || ...
            (OUT_QG_LIM == 1 && (any(gen(onld, QG) < gen(onld, QMIN) + ctol) || ...
                                any(gen(onld, QG) > gen(onld, QMAX) - ctol) || ...
                                any(gen(onld, MU_QMIN) > ptol) || ...
                                any(gen(onld, MU_QMAX) > ptol))))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Dispatchable Load Constraints                                            |');
        fprintf(fd, '\n================================================================================');
    end
    if ~isempty(onld) && (OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                             (any(gen(onld, PG) < gen(onld, PMIN) + ctol) || ...
                              any(gen(onld, PG) > gen(onld, PMAX) - ctol) || ...
                              any(gen(onld, MU_PMIN) > ptol) || ...
                              any(gen(onld, MU_PMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                  Active Power Limits');
        fprintf(fd, '\n  #     #     Pmin mu     Pmin       Pg       Pmax    Pmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_PG_LIM == 2 || (OUT_PG_LIM == 1 && ...
                        (gen(i, PG) < gen(i, PMIN) + ctol || ...
                         gen(i, PG) > gen(i, PMAX) - ctol || ...
                         gen(i, MU_PMIN) > ptol || gen(i, MU_PMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, PG) < gen(i, PMIN) + ctol || gen(i, MU_PMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, PG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [PMIN, PG, PMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [PMIN, PMAX]));
                end
                if gen(i, PG) > gen(i, PMAX) - ctol || gen(i, MU_PMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_PMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% dispatchable load Q constraints
    if ~isDC && ~isempty(onld) && (OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                             (any(gen(onld, QG) < gen(onld, QMIN) + ctol) || ...
                              any(gen(onld, QG) > gen(onld, QMAX) - ctol) || ...
                              any(gen(onld, MU_QMIN) > ptol) || ...
                              any(gen(onld, MU_QMAX) > ptol))))
        fprintf(fd, '\n Gen   Bus                 Reactive Power Limits');
        fprintf(fd, '\n  #     #     Qmin mu     Qmin       Qg       Qmax    Qmax mu');
        fprintf(fd, '\n----  -----   -------   --------  --------  --------  -------');
        for k = 1:length(onld)
            i = onld(k);
            if OUT_QG_LIM == 2 || (OUT_QG_LIM == 1 && ...
                        (gen(i, QG) < gen(i, QMIN) + ctol || ...
                         gen(i, QG) > gen(i, QMAX) - ctol || ...
                         gen(i, MU_QMIN) > ptol || gen(i, MU_QMAX) > ptol))
                fprintf(fd, '\n%4d%6d ', i, gen(i, GEN_BUS));
                if gen(i, QG) < gen(i, QMIN) + ctol || gen(i, MU_QMIN) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMIN));
                else
                    fprintf(fd, '      -   ');
                end
                if gen(i, QG)
                    fprintf(fd, '%10.2f%10.2f%10.2f', gen(i, [QMIN, QG, QMAX]));
                else
                    fprintf(fd, '%10.2f       -  %10.2f', gen(i, [QMIN, QMAX]));
                end
                if gen(i, QG) > gen(i, QMAX) - ctol || gen(i, MU_QMAX) > ptol
                    fprintf(fd, '%10.3f', gen(i, MU_QMAX));
                else
                    fprintf(fd, '      -   ');
                end
            end
        end
        fprintf(fd, '\n');
    end

    %% line flow constraints
    if upper(mpopt.opf.flow_lim(1)) == 'P' || isDC  %% |P| limit
        Ff = branch(:, PF);
        Ft = branch(:, PT);
        str = '\n  #     Bus    Pf  mu     Pf      |Pmax|      Pt      Pt  mu   Bus';
    elseif upper(mpopt.opf.flow_lim(1)) == 'I'      %% |I| limit
        Ff = abs( (branch(:, PF) + 1j * branch(:, QF)) ./ V(e2i(branch(:, F_BUS))) );
        Ft = abs( (branch(:, PT) + 1j * branch(:, QT)) ./ V(e2i(branch(:, T_BUS))) );
        str = '\n  #     Bus   |If| mu    |If|     |Imax|     |It|    |It| mu   Bus';
    else                                            %% |S| limit
        Ff = abs(branch(:, PF) + 1j * branch(:, QF));
        Ft = abs(branch(:, PT) + 1j * branch(:, QT));
        str = '\n  #     Bus   |Sf| mu    |Sf|     |Smax|     |St|    |St| mu   Bus';
    end
    if any(branch(:, RATE_A) ~= 0) && (OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
                        (any(abs(Ff) > branch(:, RATE_A) - ctol) || ...
                         any(abs(Ft) > branch(:, RATE_A) - ctol) || ...
                         any(branch(:, MU_SF) > ptol) || ...
                         any(branch(:, MU_ST) > ptol))))
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\n|     Branch Flow Constraints                                                  |');
        fprintf(fd, '\n================================================================================');
        fprintf(fd, '\nBrnch   From     "From" End        Limit       "To" End        To');
        fprintf(fd, str);
        fprintf(fd, '\n-----  -----  -------  --------  --------  --------  -------  -----');
        for i = 1:nl
            if branch(i, RATE_A) ~= 0 && (OUT_LINE_LIM == 2 || (OUT_LINE_LIM == 1 && ...
                   (abs(Ff(i)) > branch(i, RATE_A) - ctol || ...
                    abs(Ft(i)) > branch(i, RATE_A) - ctol || ...
                    branch(i, MU_SF) > ptol || branch(i, MU_ST) > ptol)))
                fprintf(fd, '\n%4d%7d', i, branch(i, F_BUS));
                if Ff(i) > branch(i, RATE_A) - ctol || branch(i, MU_SF) > ptol
                    fprintf(fd, '%10.3f', branch(i, MU_SF));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%9.2f%10.2f%10.2f', ...
                    [Ff(i), branch(i, RATE_A), Ft(i)]);
                if Ft(i) > branch(i, RATE_A) - ctol || branch(i, MU_ST) > ptol
                    fprintf(fd, '%10.3f', branch(i, MU_ST));
                else
                    fprintf(fd, '      -   ');
                end
                fprintf(fd, '%6d', branch(i, T_BUS));
            end
        end
        fprintf(fd, '\n');
    end
end

%% execute userfcn callbacks for 'printpf' stage
if have_results_struct && isfield(results, 'userfcn') && (success || OUT_FORCE)
    if ~isOPF   %% turn off option for all constraints if it isn't an OPF
        mpopt = mpoption(mpopt, 'out.lim.all', 0);
    end
    run_userfcn(results.userfcn, 'printpf', results, fd, mpopt);
end
if OUT_ANY && ~success
    if OUT_FORCE
        if isSDP
            fprintf(fd, '\n>>>>>  Solution does NOT satisfy rank and/or consistency conditions (%.2f seconds).  <<<<<\nmineigratio = %0.5g, zero_eval = %0.5g\n', et, mineigratio, zero_eval);
        else
            fprintf(fd, '\n>>>>>  Did NOT converge (%.2f seconds)  <<<<<\n', et);
        end
    end
    fprintf('\n');
end

function TorF = isload(gen)
%ISLOAD  Checks for dispatchable loads.
%   TORF = ISLOAD(GEN) returns a column vector of 1's and 0's. The 1's
%   correspond to rows of the GEN matrix which represent dispatchable loads.
%   The current test is Pmin < 0 AND Pmax == 0.
%   This may need to be revised to allow sensible specification
%   of both elastic demand and pumped storage units.

%   MATPOWER
%   Copyright (c) 2005-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

TorF = gen(:, PMIN) < 0 & gen(:, PMAX) == 0;


function [baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(casefile)
%LOADCASE   Load .m or .mat case files or data struct in MATPOWER format.
%
%   [BASEMVA, BUS, GEN, BRANCH, AREAS, GENCOST] = LOADCASE(CASEFILE)
%   [BASEMVA, BUS, GEN, BRANCH, GENCOST] = LOADCASE(CASEFILE)
%   [BASEMVA, BUS, GEN, BRANCH] = LOADCASE(CASEFILE)
%   MPC = LOADCASE(CASEFILE)
%
%   Returns the individual data matrices or a struct containing them as fields.
%
%   Here CASEFILE is either (1) a struct containing the fields baseMVA,
%   bus, gen, branch and, optionally, areas and/or gencost, or (2) a string
%   containing the name of the file. If CASEFILE contains the extension
%   '.mat' or '.m', then the explicit file is searched. If CASEFILE contains
%   no extension, then LOADCASE looks for a MAT-file first, then for an
%   M-file.  If the file does not exist or doesn't define all required
%   matrices, the routine aborts with an appropriate error message.
%
%   Alternatively, it can be called with the following syntax, though this
%   option is now deprecated and will be removed in a future version:
%
%   [BASEMVA, BUS, GEN, BRANCH, AREAS, GENCOST, INFO] = LOADCASE(CASEFILE)
%   [MPC, INFO] = LOADCASE(CASEFILE)
%
%   In this case, the function will not abort, but INFO will contain an exit
%   code as follows:
%
%       0:  all variables successfully defined
%       1:  input argument is not a string or struct
%       2:  specified extension-less file name does not exist in search path
%       3:  specified MAT-file does not exist in search path
%       4:  specified M-file does not exist in search path
%       5:  specified file fails to define all matrices or contains syntax err
%
%   If the input data is from an M-file or MAT-file defining individual
%   data matrices, or from a struct with out a 'version' field whose
%   GEN matrix has fewer than 21 columns, then it is assumed to be a
%   MATPOWER case file in version 1 format, and will be converted to
%   version 2 format.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

info = 0;
if nargout < 3
    return_as_struct = true;
else
    return_as_struct = false;
end
if nargout >= 5
    expect_gencost = true;
    if nargout > 5
        expect_areas = true;
    else 
        expect_areas = false;
    end
else
    expect_gencost = false;
    expect_areas = false;
end

%%-----  read data into struct  -----
if ischar(casefile)
    [pathstr, fname, ext] = fileparts(casefile);
    if isempty(ext)
        if exist(fullfile(pathstr, [fname '.mat']), 'file') == 2
            ext = '.mat';
        elseif exist(fullfile(pathstr, [fname '.m']), 'file') == 2
            ext = '.m';
        else
            info = 2;
        end
    end
    
    %% attempt to read file
    if info == 0
        if strcmp(ext,'.mat')       %% from MAT file
            try
                s = load(fullfile(pathstr, fname));
                if isfield(s, 'mpc')    %% it's a struct
                    s = s.mpc;
                else                    %% individual data matrices
                    s.version = '1';
                end
            catch
                info = 3;
            end
        elseif strcmp(ext,'.m')     %% from M file
            if ~isempty(pathstr)
                cwd = pwd;          %% save working directory to string
                cd(pathstr);        %% cd to specified directory
            end
            try                             %% assume it returns a struct
                s = feval(fname);
            catch
                info = 4;
            end
            if info == 0 && ~isstruct(s)    %% if not try individual data matrices
                clear s;
                s.version = '1';
                if expect_gencost
                    try
                        [s.baseMVA, s.bus, s.gen, s.branch, ...
                            s.areas, s.gencost] = feval(fname);
                    catch
                        info = 4;
                    end
                else
                    if return_as_struct
                        try
                            [s.baseMVA, s.bus, s.gen, s.branch, ...
                                s.areas, s.gencost] = feval(fname);
                        catch
                            try
                                [s.baseMVA, s.bus, s.gen, s.branch] = feval(fname);
                            catch
                                info = 4;
                            end
                        end
                    else
                        try
                            [s.baseMVA, s.bus, s.gen, s.branch] = feval(fname);
                        catch
                            info = 4;
                        end
                    end
                end
            end
            if info == 4 && exist(fullfile(pathstr, [fname '.m']), 'file') == 2
                info = 5;
                err5 = lasterr;
            end
            if ~isempty(pathstr)    %% change working directory back to original
                cd(cwd);
            end
        end
    end
elseif isstruct(casefile)
    s = casefile;
else
    info = 1;
end

%%-----  check contents of struct  -----
if info == 0
    %% check for required fields
    if expect_areas && ~isfield(s,'areas')
        s.areas = [];   %% add empty missing areas if needed for output
    end
    if ~( isfield(s,'baseMVA') && isfield(s,'bus') && ...
            isfield(s,'gen') && isfield(s,'branch') ) || ...
            ( expect_gencost && ~isfield(s, 'gencost') )
        info = 5;           %% missing some expected fields
        err5 = 'missing data';
    else
        %% remove empty areas if not needed
        if isfield(s, 'areas') && isempty(s.areas) && ~expect_areas
            s = rmfield(s, 'areas');
        end

        %% all fields present, copy to mpc
        mpc = s;
        if ~isfield(mpc, 'version') %% hmm, struct with no 'version' field
            if size(mpc.gen, 2) < 21    %% version 2 has 21 or 25 cols
                mpc.version = '1';
            else
                mpc.version = '2';
            end
        end
        if strcmp(mpc.version, '1')
            % convert from version 1 to version 2
            [mpc.gen, mpc.branch] = mpc_1to2(mpc.gen, mpc.branch);
            mpc.version = '2';
        end
    end
end

%%-----  define output variables  -----
if return_as_struct
    bus = info;
end

if info == 0    %% no errors
    if return_as_struct
        baseMVA = mpc;
    else
        baseMVA = mpc.baseMVA;
        bus     = mpc.bus;
        gen     = mpc.gen;
        branch  = mpc.branch;
        if expect_gencost
            if expect_areas
                areas   = mpc.areas;
                gencost = mpc.gencost;
            else
                areas = mpc.gencost;
            end
        end
    end
else            %% we have a problem captain
    if nargout == 2 || nargout == 7   %% return error code
        if return_as_struct
            baseMVA = struct([]);
        else
            baseMVA = []; bus = []; gen = []; branch = [];
            areas = []; gencost = [];
        end
    else                                            %% die on error
        switch info
            case 1,
                error('loadcase: input arg should be a struct or a string containing a filename');
            case 2,
                error('loadcase: specified case not in MATLAB''s search path');
            case 3,
                error('loadcase: specified MAT file does not exist');
            case 4,
                error('loadcase: specified M file does not exist');
            case 5,
                error('loadcase: syntax error or undefined data matrix(ices) in the file\n%s', err5);
            otherwise,
                error('loadcase: unknown error');
        end
    end
end



function [gen, branch] = mpc_1to2(gen, branch)

%% define named indices into bus, gen, branch matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  gen  -----
%% use the version 1 values for column names
if size(gen, 2) > APF
    error('mpc_1to2: gen matrix appears to already be in version 2 format');
end
shift = MU_PMAX - PMIN - 1;
tmp = num2cell([MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] - shift);
[MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = deal(tmp{:});

%% add extra columns to gen
tmp = zeros(size(gen, 1), shift);
if size(gen, 2) >= MU_QMIN
    gen = [ gen(:, 1:PMIN) tmp gen(:, MU_PMAX:MU_QMIN) ];
else
    gen = [ gen(:, 1:PMIN) tmp ];
end

%%-----  branch  -----
%% use the version 1 values for column names
shift = PF - BR_STATUS - 1;
tmp = num2cell([PF, QF, PT, QT, MU_SF, MU_ST] - shift);
[PF, QF, PT, QT, MU_SF, MU_ST] = deal(tmp{:});

%% add extra columns to branch
tmp = ones(size(branch, 1), 1) * [-360 360];
tmp2 = zeros(size(branch, 1), 2);
if size(branch, 2) >= MU_ST
    branch = [ branch(:, 1:BR_STATUS) tmp branch(:, PF:MU_ST) tmp2 ];
elseif size(branch, 2) >= QT
    branch = [ branch(:, 1:BR_STATUS) tmp branch(:, PF:QT) ];
else
    branch = [ branch(:, 1:BR_STATUS) tmp ];
end
