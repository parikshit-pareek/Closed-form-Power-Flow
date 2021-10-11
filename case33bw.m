function mpc = case33bw
%CASE33BW  Power flow data for 33 bus distribution system from Baran & Wu
%    Please see CASEFORMAT for details on the case file format.
%
%    Data from ...
%       M. E. Baran and F. F. Wu, "Network reconfiguration in distribution
%       systems for loss reduction and load balancing," in IEEE Transactions
%       on Power Delivery, vol. 4, no. 2, pp. 1401-1407, Apr 1989.
%       doi: 10.1109/61.25627
%       URL: http://doi.org/10.1109/61.25627

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [  %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	12.66	1	1	1;
	2	1	100	60	0	0	1	1	0	12.66	1	1.1	0.9;
	3	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	4	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	5	1	60	30	0	0	1	1	0	12.66	1	1.1	0.9;
	6	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	7	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	8	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	9	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	10	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	11	1	45	30	0	0	1	1	0	12.66	1	1.1	0.9;
	12	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	13	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	14	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	15	1	60	10	0	0	1	1	0	12.66	1	1.1	0.9;
	16	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	17	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	18	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	19	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	20	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	21	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	22	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	23	1	90	50	0	0	1	1	0	12.66	1	1.1	0.9;
	24	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	25	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	26	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	27	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	28	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	29	1	120	70	0	0	1	1	0	12.66	1	1.1	0.9;
	30	1	200	600	0	0	1	1	0	12.66	1	1.1	0.9;
	31	1	150	70	0	0	1	1	0	12.66	1	1.1	0.9;
	32	1	210	100	0	0	1	1	0	12.66	1	1.1	0.9;
	33	1	60	40	0	0	1	1	0	12.66	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1.0	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.0922	0.0470	0	0	0	0	0	0	1	-360	360;
	2	3	0.4930	0.2511	0	0	0	0	0	0	1	-360	360;
	3	4	0.3660	0.1864	0	0	0	0	0	0	1	-360	360;
	4	5	0.3811	0.1941	0	0	0	0	0	0	1	-360	360;
	5	6	0.8190	0.7070	0	0	0	0	0	0	1	-360	360;
	6	7	0.1872	0.6188	0	0	0	0	0	0	1	-360	360;
	7	8	0.7114	0.2351	0	0	0	0	0	0	1	-360	360;
	8	9	1.0300	0.7400	0	0	0	0	0	0	1	-360	360;
	9	10	1.0440	0.7400	0	0	0	0	0	0	1	-360	360;
	10	11	0.1966	0.0650	0	0	0	0	0	0	1	-360	360;
	11	12	0.3744	0.1238	0	0	0	0	0	0	1	-360	360;
	12	13	1.4680	1.1550	0	0	0	0	0	0	1	-360	360;
	13	14	0.5416	0.7129	0	0	0	0	0	0	1	-360	360;
	14	15	0.5910	0.5260	0	0	0	0	0	0	1	-360	360;
	15	16	0.7463	0.5450	0	0	0	0	0	0	1	-360	360;
	16	17	1.2890	1.7210	0	0	0	0	0	0	1	-360	360;
	17	18	0.7320	0.5740	0	0	0	0	0	0	1	-360	360;
	2	19	0.1640	0.1565	0	0	0	0	0	0	1	-360	360;
	19	20	1.5042	1.3554	0	0	0	0	0	0	1	-360	360;
	20	21	0.4095	0.4784	0	0	0	0	0	0	1	-360	360;
	21	22	0.7089	0.9373	0	0	0	0	0	0	1	-360	360;
	3	23	0.4512	0.3083	0	0	0	0	0	0	1	-360	360;
	23	24	0.8980	0.7091	0	0	0	0	0	0	1	-360	360;
	24	25	0.8960	0.7011	0	0	0	0	0	0	1	-360	360;
	6	26	0.2030	0.1034	0	0	0	0	0	0	1	-360	360;
	26	27	0.2842	0.1447	0	0	0	0	0	0	1	-360	360;
	27	28	1.0590	0.9337	0	0	0	0	0	0	1	-360	360;
	28	29	0.8042	0.7006	0	0	0	0	0	0	1	-360	360;
	29	30	0.5075	0.2585	0	0	0	0	0	0	1	-360	360;
	30	31	0.9744	0.9630	0	0	0	0	0	0	1	-360	360;
	31	32	0.3105	0.3619	0	0	0	0	0	0	1	-360	360;
	32	33	0.3410	0.5302	0	0	0	0	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];


%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
Vbase = mpc.bus(1, BASE_KV) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);

%% convert loads from kW to MW
mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;



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

