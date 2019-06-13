% 20180701
% Generate OTOC using OperatorClass
% y= Tr(|rho(t),O|^2)
% arguments: rho, H, O and tlist
% rho: the time evolved operator
% H: Hamiltonian
% O: the static operator
% tlist: time to evaluate OTOC
% created by Pai Peng
function y=OTOC(rho, H, O, tlist)
if ~IsSameSym({rho,H,O})
    error('The system size or symmetry class of the input operators are not the same.')
end
if isempty(H.eigsys)
    H.diagonalize
end

rhoV=cell(1,length(rho.matrix));
OV=cell(1,length(rho.matrix));

for p=1:length(rho.matrix)
    rhoV{p}=H.eigsys{p}.V'*rho.matrix{p}*H.eigsys{p}.V;
    OV{p}=H.eigsys{p}.V'*O.matrix{p}*H.eigsys{p}.V;
end

y=zeros(length(tlist),length(rho.matrix));
for q=1:length(tlist)
    t=tlist(q);
    for p=1:length(rho.matrix)
        Ut=sparse(diag(exp(-1i*H.eigsys{p}.D*t)));
        rhoVt=Ut*rhoV{p}*Ut';
        y(q,p)=sum(sum(abs(rhoVt*OV{p}-OV{p}*rhoVt).^2));
    end
end

y=sum(y,2);