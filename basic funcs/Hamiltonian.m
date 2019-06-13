% 20180423
% Hamiltonian generation
% input:
% # of atoms
% bc: boundary condition, 'o' for obc, 'p' for pbc
% cplist: gives coupling range and strength
% type: 'dipolar' or 'Ising' or 'Heisenbergs' or a 3-element vector
% specifying the prefactor of XX, YY and ZZ.
% direction: 1 for x, etc
% Created by Pai Peng
% 20180814
% Includes double quantum Hamiltonian XX-YY with type code 'dq'
% by PP
function H_int=Hamiltonian(N_atom,bc,cplist,type,direction)
switch nargin % check the validity of # of arguments and direction
    case 4
        if ~strcmp(type,'Heisenberg') && ~strcmp(type,'dq') && ~isa(type,'double')
            error('Specify direction')
        end
    case 5
        if ~ismember(direction,[1,2,3])
            error('Invalid direction')
        end
        
    otherwise
        error('Invalid # of inputs')
end

H_int=sparse(2^N_atom,2^N_atom);

for pp=1:length(cplist)
    for p=1:N_atom
        switch bc
            case 'p'
                q=mod(p+pp-1,N_atom)+1;
            case 'o'
                q=p+pp;
                if q>N_atom
                    continue
                end
            otherwise
                error('Invalid boundary condition')
        end
        if isa(type,'double')
            H_int=H_int+...
                cplist(pp)*(type(1)*Kron2body(N_atom,1,1,p,q)...
                +type(2)*Kron2body(N_atom,2,2,p,q)...
                +type(3)*Kron2body(N_atom,3,3,p,q));
        else
            switch type
                case 'dipolar'
                    H_int=H_int+0.5*cplist(pp)*...
                        (3*Kron2body(N_atom,direction,direction,p,q)-...
                        (Kron2body(N_atom,1,1,p,q)+Kron2body(N_atom,2,2,p,q)+Kron2body(N_atom,3,3,p,q)));
                case 'Ising'
                    H_int=H_int+cplist(pp)*Kron2body(N_atom,direction,direction,p,q);
                case 'Heisenberg'
                    H_int=H_int+...
                        cplist(pp)*(Kron2body(N_atom,1,1,p,q)+Kron2body(N_atom,2,2,p,q)+Kron2body(N_atom,3,3,p,q));
                case 'dq'
                    H_int=H_int+...
                        cplist(pp)*(Kron2body(N_atom,1,1,p,q)-Kron2body(N_atom,2,2,p,q));
                otherwise
                    error('Invalid interaction type')
            end
        end
    end
end