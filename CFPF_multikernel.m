%% Code CFPF Learning using Different Kernels 

% Paper: 

% Date of Code: 26 Apr 2020, Covid-19 Lockdown time, 4:38 in the morning
% Copyright: Parikshit Pareek, NTU Singapore, pare0001@ntu.edu.sg

% % % ------------------------------------------------------------------------------------------------
clear
clc


data = case33bw; % The system under consideration


N_training = 60; % number of samples in the tarining lattice structure 
N_testing = 100; % number of samples for testing


scale=1;
pqbus = data.bus(data.bus(:,2)==1,1); %% PQ Buses where voltage is unknown
pvbus = data.bus(data.bus(:,2)==2,1); %% PV Buses
pvqbus = data.bus(data.bus(:,2)<3,1); %% PV and PQ bus
lbus = data.bus(data.bus(:,3)> 0,1);
nbus = length(data.bus(:,1));
%% --------------------------- Type of uncertainty ----------------------------------------------
%  
% All bus with nonzero load percentage Uncertainty 
% # : >  PV+PQ buses have uncertain Injection of P and  
% # : > PQ buses have uncertain Injection of Q
% Uncertainty level
sfP=0.2; % Real power percentage fracation of base-load
sfQ=sfP; % Reactive power percentage fracation of base-load

parfor kr=1:4

[D_learningkr{kr},Loutkr{kr},Base]=IPF_learning_polykernel(data,pvbus,pqbus,N_training,N_testing,kr,sfP,sfQ,scale);
D_learning=D_learningkr{kr};
Lout = Loutkr{kr}; 

[~,Mcskr{kr}]=MCS_output(data,D_learning.xs,D_learning.rbus,D_learning.D,N_testing,Lout,scale,pqbus);
Mcs = Mcskr{kr};

%% Error Analysis
% Voltage Magnitude
Vac=Mcs.V;
erV=abs(Vac-Lout.muV);
maeV(:,kr) = sum(abs(erV))/N_testing;

% Voltage Angle
Thac=Mcs.Thac;
erTh=abs(Mcs.Thac-Lout.muTh);
maeTh(:,kr) = sum(abs(erTh))/N_testing;



end


%% Function for learning at different kernel's 
% For learning all PQ+PV buses theta and PQ bus |V| 

% Date 20 Apr 2020: Modified on 26 Apr
% Parikshit Pareek @NTUSg


function [D_learning,Lout,Base]=IPF_learning_polykernel(data,~,pqbus,N,nt,kr,sfP,sfQ,scale)
% nc:= Number of control variables
% xclimit:= Control limits for all the D dimensions
% scale := Base for voltage learning
% ==============       Output     =============================================================
%       Lout := All output related to the learning of voltage maginitude function
% D_learning := Learning input dataset
%       Base := Base case results
%       Bopt := Optimal Bound Results
% ResControl := Voltage Control Results
pvqbus = data.bus(data.bus(:,2)<3,1);
% nbus = length(data.bus(:,1));

% pqbus = data.bus(data.bus(:,2)==1,1); %% PQ Buses where voltage is unknown

%%-     Construction of input set load uncertainty input dimension
%   -- Only load uncertainty without the powerfactor constraint  -----------
% [data_up,xx,xs,xlimit,D,rbus]=input_dataset_NetInjection(data,Tu,N,nt,2,sfP,sfQ,pqbus,pvbus);
[data_up,xx,xs,xlimit,D,rbus]=input_dataset_Load(data,1,N,nt,2,sfP,sfQ,pqbus);



%------------------ Generating Target for Training Input xx ----------------------------------------
% D_learning = Learning Dataset
parfor i=1:N
[V(:,:,i),Sg(:,:,i),Sij(:,:,i),sol{i},J{i},Y(:,:,i)]=Sampling_Jaco(data_up{i});
end

D_learning.V=V; D_learning.Sg=Sg; D_learning.Sij=Sij; 
D_learning.xx=xx; D_learning.xs=xs; D_learning.xlimit=xlimit; D_learning.D=D; D_learning.sol=sol; 
D_learning.J=J;  D_learning.data_up=data_up; D_learning.rbus=rbus;
D_learning.Y=Y;

%% Specify the mean, covariance and likelihood functions 
meanfunc = [];                                                  % meanZero
likfunc = @likGauss;                                            % Gaussian likelihood 
if kr == 1
    covfunc = {@covPoly,1}; % Linear (Polynomial of degree 1)  
elseif kr == 2
    covfunc = {@covPoly,2}; % Quadratic kernel (Polynomial of degree 2)  
elseif kr == 3
        covfunc = @covRQiso;  % Squared Exponental covariance function @covSEiso
elseif kr == 4
        covfunc = @covSEiso;  % Squared Exponental covariance function @covSEiso
elseif kr == 5
    covfunc = {@covMaterniso,3}; % Polynomial of degree 5
end
%% ======= Using parallel programming for learning Voltage magnitude ==========
parfor k=1:length(pqbus)
    yy=reshape(V(pqbus(k),1,:),[N,1])*scale;
    ytV(:,k)=yy;
        if kr < 4
            hyp = struct('mean', [], 'cov', [1; -1;1], 'lik', -1);
        else
            hyp = struct('mean',[], 'cov', [5; -1], 'lik', -1);
        end
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xx, yy);
%     hyp_store{k}=hyp;
    sf_lV(:,k)=exp(hyp.cov);
    [muV(:,k),s2V(:,k),~,~,~,post] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy,xs);
    alphaV(:,k)=post.alpha;
end

% Defining a structure Lopt: Learning Output 
Lout.alphaV =alphaV;  Lout.muV = muV; % Lout.Smax_V=max_V;  Lout.Smin_V=min_V;  
  Lout.s2V = s2V; Lout.sf_lV=sf_lV;  Lout.ytV = ytV; %Lout.K=K; Lout.btaV=btaV; Lout.rkhs_sqV=rkhs_sqV; 

% Standarad Form of Linear Function of voltage: V= W*s + C
if kr==1
      for k=1:length(pqbus)
                % Making the constraint 
                bV(:,:,k)=xx.*repmat(alphaV(:,k),[1,D]);
                wV(:,k)= sum(bV(:,:,k));
                WV(:,k) = (sf_lV(3,k)*sf_lV(3,k)*wV(:,k))/(sf_lV(1,k)^2);
                CV(k)=sf_lV(3,k)*sf_lV(3,k)*sf_lV(2,k)*sum(alphaV(:,k));                          
      end 
      WV=WV'; CV=CV';

Lout.WV=WV;
Lout.CV=CV;
end

%% ======= Using parallel programming for learning Voltage solution ==========
parfor k=1:length(pvqbus)
    yy=reshape(V(pvqbus(k),2,:),[N,1]);
    ytTh(:,k)=yy;
        if kr < 4
            hyp = struct('mean', [], 'cov', [1; -1;1], 'lik', -1);
        else
            hyp = struct('mean',[], 'cov', [5; -1], 'lik', -1);
        end
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xx, yy);
%     hyp_store{k}=hyp;
    sf_lTh(:,k)=exp(hyp.cov);
    [muTh(:,k),s2Th(:,k),~,~,~,post] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy,xs);

    alphaTh(:,k)=post.alpha;
end

% Updating a structure Lopt: Learning Output 
Lout.alphaTh =alphaTh; Lout.muTh = muTh; %Lout.btaTh=btaTh;   Lout.Smax_Th=max_Th;  Lout.Smin_Th=min_Th;  
  Lout.s2Th = s2Th; Lout.sf_lTh=sf_lTh;  Lout.ytTh = ytTh; %Lout.Kt=Kt; Lout.rkhs_sqTh=rkhs_sqTh; 



[Base.V,Base.Sg,Base.Sij,Base.sol,Base.J]=Sampling_Jaco(data);

end

