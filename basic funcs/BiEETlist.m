% 20180812
% Generate bipartite entanglement entropy using OperatorClass
% arguments: rho, H, O and tlist
% rho: the time evolved operator
% H: Hamiltonian
% tlist: time to evaluate OTOC
% created by Pai Peng
function y=BiEETlist(rho, H, tlist, n)
if ~IsSameSym({rho,H})
    error('The system size or symmetry class of the input operators are not the same.')
end
if ~isempty(rho.sym)
    error('Operator must not be symmetrized')
end
if isempty(H.eigsys)
    H.diagonalize
end

rhoV=H.eigsys{1}.V'*rho.matrix{1}*H.eigsys{1}.V;

y=zeros(length(tlist),1);
for q=1:length(tlist)
    t=tlist(q);
    Ut=sparse(diag(exp(-1i*H.eigsys{1}.D*t)));
    rhot=H.eigsys{1}.V*Ut*rhoV*Ut'*H.eigsys{1}.V';
    rhott=OperatorClass(rho.L);% OperatorClass of rhot
    rhott.matrix={rhot};
    y(q)=BipartiteEE(rhott,n);
end
