%20170912
%returns individual pauli matrices of a N-atom array
%1 argument N: returns all 3N paulis
%3 arguments N direction position: returns the Pauli for specified
%direction and position
function y=IndivPauliSparse(N,direction,position)
switch nargin
    case 1
        Sigma_xi=cell(1,N);
        Sigma_yi=cell(1,N);
        Sigma_zi=cell(1,N);
        
        Sigma_xi{1}=sparse(sigma_x);
        Sigma_yi{1}=sparse(sigma_y);
        Sigma_zi{1}=sparse(sigma_z);
        for p=2:N
            Sigma_xi{p}=sparse(eye(2));
            Sigma_yi{p}=sparse(eye(2));
            Sigma_zi{p}=sparse(eye(2));
        end
        
        for p=2:N
            for pp=1:N
                if p==pp
                    Sigma_xi{pp}=kron(Sigma_xi{pp},sparse(sigma_x));
                    Sigma_yi{pp}=kron(Sigma_yi{pp},sparse(sigma_y));
                    Sigma_zi{pp}=kron(Sigma_zi{pp},sparse(sigma_z));
                else
                    Sigma_xi{pp}=kron(Sigma_xi{pp},sparse(eye(2)));
                    Sigma_yi{pp}=kron(Sigma_yi{pp},sparse(eye(2)));
                    Sigma_zi{pp}=kron(Sigma_zi{pp},sparse(eye(2)));
                end
            end
        end
        
        y={Sigma_xi,Sigma_yi,Sigma_zi};
    case 3
        y=sparse(1);
        for p=1:position-1
            y=kron(y,sparse(eye(2)));
        end
        switch direction
            case 1
                y=kron(y,sparse(sigma_x));
            case 2
                y=kron(y,sparse(sigma_y));
            case 3
                y=kron(y,sparse(sigma_z));
            otherwise
                error('Invalid direction')
        end
        for p=position+1:N_atom
            y=kron(y,sparse(eye(2)));
        end
    otherwise
        error('Invalid # of inputs')
end
