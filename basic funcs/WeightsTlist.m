% 20181005
% Generate correlation length for tlist using OperatorClass
% arguments: rho, H, and tlist
% rho: the time evolved operator
% H: Hamiltonian
% tlist: time to evaluate correlation length
% qu: cell, containing 'CL' and/or 'Hamming' and/or number m (for notm)
% created by Pai Peng
function y=WeightsTlist(rho, H, tlist, qu)
if ~IsSameSym({rho,H})
    error('The system size or symmetry class of the input operators are not the same.')
end
if ~isempty(rho.sym)
    error('rho must not be symmetrized.')
end
if isempty(H.eigsys)
    H.diagonalize
end

rhoV=cell(1,length(rho.matrix));


    rhoV{1}=H.eigsys{1}.V'*rho.matrix{1}*H.eigsys{1}.V;


y=zeros(length(tlist),length(qu));
for q=1:length(tlist)
    t=tlist(q);
    for p=1
        Ut=sparse(diag(exp(-1i*H.eigsys{p}.D*t)));
        rhot=H.eigsys{p}.V*Ut*rhoV{p}*Ut'*H.eigsys{p}.V';
    end
    Cl=Weights(rhot,'sum','square',qu);
    y(q,:)=sum(Cl.*(0:(size(Cl,1)-1)).',1);
end
