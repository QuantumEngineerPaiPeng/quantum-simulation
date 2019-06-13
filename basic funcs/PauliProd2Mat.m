% 20180731
% Convert a Pauli product operator to a matrix
% Input: e.g. '11eg+-xyz'
% Output: 2^L by 2^L matrix
% Created by Pai Peng
function y=PauliProd2Mat(x)
L=length(x);
y=Pauli(x(1));
for p=2:L
    y=kron(y,Pauli(x(p)));
end