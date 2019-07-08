% 20180630
% Generate local Pauli operators in OperatorClass
% Input: 
% L: system size
% n: Pauli on nth site being
% tag: 'x', 'y' 'z' '+' or '-' giving the Pauli operator
% Created by Pai Peng
function y=LocalPauli(L,n,tag)
switch tag
    case 'x'
        pauli=sparse(sigma_x);
    case 'y'
        pauli=sparse(sigma_y);
    case 'z'
        pauli=sparse(sigma_z);
    case '-'
        pauli=sparse([0,1;0,0]);
    case '+'
        pauli=sparse([0,0;1,0]);
    otherwise
        error('3rd input must be ''x'', ''y'' ''z'' ''+'' or ''-''')
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