% 20180723
% for a given Hamiltonian, get the canonical ensemble averaged energy
% input: inverse temperature beta, Hamiltonian H
% output: y:energy
% created by Pai Peng
function y=CanEnsEnergy(beta,H)
if isempty(H.eigsys)
   H.diagonalize()
end
E=0;
Z=0;
for p=1:length(H.matrix)
    E=E+H.eigsys{p}.D.'*exp(-beta*H.eigsys{p}.D);
    Z=Z+sum(exp(-beta*H.eigsys{p}.D));
end

y=E/Z;