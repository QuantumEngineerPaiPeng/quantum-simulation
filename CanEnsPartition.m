% 20180723
% for a given Hamiltonian, get the canonical ensemble partition function
% input: inverse temperature beta, Hamiltonian H
% output: Z:partition function
% created by Pai Peng
function Z=CanEnsPartition(beta,H)
if isempty(H.eigsys)
   H.diagonalize()
end
Z=0;
for p=1:length(H.matrix)
    Z=Z+sum(exp(-beta*H.eigsys{p}.D));
end