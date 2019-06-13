% 20180920
% function to fake a MBL like Hamiltonian in tau_z basis
% Input: L: system size, gamma: coefficient decay rate
% Created by Pai Peng
function y=fakeMBLHam(L,gamma)
CL=CLcell(L,2);
y=zeros(L+1,1);
for p=3:length(CL)
%     y(CL{p})=exp(normrnd(-gamma*(p-3),1.5*log(10),length(CL{p}),1));
%     y(CL{p})=exp(normrnd(-gamma*(p-3),0.3*gamma*(p-3),length(CL{p}),1));
%     y(CL{p})=exp(normrnd(-2*log(10),0.8*log(10),length(CL{p}),1));
    y(CL{p})=exp(-gamma*(p-1));
%     y(CL{p})=y(CL{p}).*(1+1e-1*rand(size(y(CL{p}))));
    if p>3
        y(CL{p})=y(CL{p}).*(sign(rand(size(y(CL{p})))-0.5));
    end
%     y(CL{p})=y(CL{p}).*(rand(size(y(CL{p})))-0.5);
end
y(CL{2})=10*(rand(size(y(CL{2}))))*2;
y=reconZ(y);