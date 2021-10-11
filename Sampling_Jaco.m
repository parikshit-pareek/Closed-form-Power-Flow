% Cover to sample power flow results using Matpowe 'runpf' code
% Parikshit @NTUSg


function [V,Sg,Sij,sol,J,Ybus]=Sampling_Jaco(data)

[~, bus, gen, branch, success, ~,J,Ybus] = ...
                runpf_complete(data);
Ybus=full(Ybus);
if success==1            
    V(:,1:2)  = bus(:,8:9); % V = [Vmag Theta]
    Sg(:,1:2) = gen(:,2:3); % Sg = [Pg Qg]
    Sij(:,1:2)= branch(:,1:2); % From node To node
    Sij(:,3:6)= branch(:,14:17); % Sij=[Fbus Tbus Pij Qij Pji Qji]
    t = Sij(:,3)+Sij(:,5)+1j*(Sij(:,4)+Sij(:,6)); % Apperent Power Loss
    Sij(:,7)= abs(t);
    Sij(:,8)= sqrt(Sij(:,3).*Sij(:,3)+Sij(:,4).*Sij(:,4)); % Apperent Power flow
    Sij(:,9)= sqrt(Sij(:,5).*Sij(:,5)+Sij(:,6).*Sij(:,6)); % Apperent Power flow
    Sij(:,10) = max(Sij(:,8),Sij(:,9));
    sol.bus=bus; sol.branch=branch; sol.gen=gen; sol.success=success; % structure of results
else
    V=[];Sg=[];Sij=[];sol=[];
end
sol.bus=bus; sol.branch=branch; sol.gen=gen; sol.success=success;
            
            
end