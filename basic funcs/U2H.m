% 20186021
% Generate Floquet Hamiltonian using an unitary operator
% Output: a cell containing the Hamiltonian in each symmetry sector
% Input: 
% Unitary, OperatorClass object
% Time, scalar
% created by Pai Peng
function y=U2H(U,t)
if ~isscalar(t)
    error('Time must be given as a scalar')
end
if isempty(U.eigsys)
    U.diagonalize
end
y=copy(U);
for p=1:length(U.matrix)
    y.eigsys{p}.D=1i*log(U.eigsys{p}.D)/t;
    y.matrix{p}=U.eigsys{p}.V*sparse(diag(y.eigsys{p}.D))*U.eigsys{p}.V';
end
