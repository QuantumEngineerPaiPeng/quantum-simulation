% 20180630
% Generate local Pauli operators in OperatorClass
% Input: 
% L: system size
% n: Pauli on nth site being
% tag: 'x', 'y' or 'z' giving the Pauli operator
% Created by Pai Peng
function y=LocalPauli(L,n,tag)
switch tag
    case 'x'
        pauli=sparse(sigma_x);
    case 'y'
        pauli=sparse(sigma_y);
    case 'z'
        pauli=sparse(sigma_z);
    otherwise
        error('3rd input must be ''x'', ''y'' or ''z''')
end

y=OperatorClass(L);
matr=sparse(1);
for p=1:n-1
    matr=kron(matr,speye(2));
end
matr=kron(matr,pauli);
for p=n+1:L
    matr=kron(matr,speye(2));
end
y.matrix={matr};