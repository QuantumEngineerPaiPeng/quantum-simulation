% 20180629
% Generate random field given by an operator class
% Input:
% L: system size
% randM: 3xL matrix, giving the field strength at each site along each
% direction
function y=randOpe(L,randM)
y=OperatorClass(L);
spx=sparse(sigma_x);
spy=sparse(sigma_y);
spz=sparse(sigma_z);
matr=randM(1,1)*spx+randM(2,1)*spy+randM(3,1)*spz;
for p=2:L
    matr=kron(matr,speye(2))+kron(speye(2^(p-1)),randM(1,p)*spx+randM(2,p)*spy+randM(3,p)*spz);
end
y.matrix={sparse(matr)};