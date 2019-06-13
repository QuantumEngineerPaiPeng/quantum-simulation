%20170912
%returns individual pauli matrices of a N-atom array
function y=IndivPauli(N)
Sigma_xi=cell(1,N);
Sigma_yi=cell(1,N);
Sigma_zi=cell(1,N);

Sigma_xi{1}=sigma_x;
Sigma_yi{1}=sigma_y;
Sigma_zi{1}=sigma_z;
for p=2:N
    Sigma_xi{p}=eye(2);
    Sigma_yi{p}=eye(2);
    Sigma_zi{p}=eye(2);
end

for p=2:N
    for pp=1:N
        if p==pp
            Sigma_xi{pp}=kron(Sigma_xi{pp},sigma_x);
            Sigma_yi{pp}=kron(Sigma_yi{pp},sigma_y);
            Sigma_zi{pp}=kron(Sigma_zi{pp},sigma_z);
        else
            Sigma_xi{pp}=kron(Sigma_xi{pp},eye(2));
            Sigma_yi{pp}=kron(Sigma_yi{pp},eye(2));
            Sigma_zi{pp}=kron(Sigma_zi{pp},eye(2));
        end
    end
end

y={Sigma_xi,Sigma_yi,Sigma_zi};
end
            