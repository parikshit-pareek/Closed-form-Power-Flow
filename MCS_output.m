%% MCS Function for all the buses
% Date 20 Apr 2020

%% MCS Output and error analysis for |V| at PQ buses and Theta at PV and PQ bus
function [data_upxs,Mcs]=MCS_output(data,xs,rbus,D,nt,Lout,scale,pqbus)
pvqbus = data.bus(data.bus(:,2)<3,1);
data_upxs=cell(1,nt);
scaleth=1;
% nbus = length(data.bus(:,1));

for j=1:nt
        data.bus(rbus,3:4)= [xs(j,1:D/2)' xs(j,D/2+1:end)'];
        data_upxs{j}=data;
end

 parfor i=1:nt
      [Vmc(:,:,i),Sgmc{i},Sijmc{i},solmc{i}]=Sampling_Jaco(data_upxs{i});
 end
 
 %% Error Analysis 
parfor k=1:length(pqbus)
 Vac(:,k)=Vmc(pqbus(k),1,:);
end
parfor k=1:length(pvqbus)
     Thac(:,k)=Vmc(pvqbus(k),2,:);
end
%% Error analysis of V
 Mcs.V=Vac*scale;
 Vac=Vac*scale;
 Mcs.erV_par=(abs(Vac-Lout.muV)./Vac)*100;
 
Mcs.erV_L1 = (norm(Vac*scale-Lout.muV,1)/norm(Vac*scale,1))*100;
Mcs.erV_L2 = (norm(Vac*scale-Lout.muV,2)/norm(Vac*scale,2))*100;
Mcs.erV_Linf = (norm(Vac*scale-Lout.muV,inf)/norm(Vac*scale,inf))*100;



%% Error in Theta 
Mcs.Thac=Thac*scaleth;
Thac=Thac*scaleth;
Mcs.erTh_par=((abs(Thac)*scaleth-abs(Lout.muTh))./(abs(Thac)))*100;
 
Mcs.erTh_L1 = (norm(abs(Thac)*scaleth-abs(Lout.muTh),1)/norm(abs(Thac)*scaleth,1))*100;
Mcs.erTh_L2 = (norm(abs(Thac)*scaleth-abs(Lout.muTh),2)/norm(abs(Thac)*scaleth,2))*100;
Mcs.erTh_Linf = (norm(abs(Thac)*scaleth-abs(Lout.muTh),inf)/norm(abs(Thac)*scaleth,inf))*100;




end







