% 20186021
% Generate unitary using a Hamiltonain operator
% Output: a cell containing the unitary in each symmetry sector
% Input: 
% Hamiltonian, OperatorClass object
% Time, scalar
% created by Pai Peng
function y=H2U(H,t)
if ~isscalar(t)
    error('Time must be given as a scalar')
end
if isempty(H.eigsys)
    H.diagonalize
end
y=copy(H);
% y.L=9;
for p=1:length(H.matrix)
    y.eigsys{p}.D=exp(-1i*H.eigsys{p}.D*t);
    y.matrix{p}=H.eigsys{p}.V*sparse(diag(y.eigsys{p}.D))*H.eigsys{p}.V';
end
