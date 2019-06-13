function psi0 = initial_psi( ...,
    psi_type, N_atom, psi_N ...,
)
% 20190313
    % 'up': all up, 'up/down': each spin either up or down, product state,
     % 'rand': total random in 2^N Hilbert space, 
     %'rand_pro': each spin rand direction, different spin no entangle
% Chao Yin

%% test state of IPR
if strcmp( psi_type, 'rand')
    psi0 = 2*rand(2^N_atom ,psi_N)-1 + i* ( 2*rand(2^N_atom ,psi_N)-1 ); 
    for psi_index = 1:psi_N
        psi0(:, psi_index) = psi0(:, psi_index) / norm( psi0(:, psi_index) );
    end
elseif strcmp( psi_type, 'up')
    psi0 = zeros(2^N_atom , psi_N);
    psi0(1, :) = 1;
    
elseif strcmp( psi_type, 'up/down')
    psi0 = zeros(2^N_atom , psi_N);
    for psi_index = 1: psi_N % N+1- psi_index spins: up, other psi_index-1 : down
        temp = 1;
        this_perm = randperm(N_atom);
        for N_index = 1:N_atom
            if this_perm(N_index) >= psi_index
                temp = kron(temp, [1;0]);
            else
                temp = kron(temp, [0;1]);
            end
        end
        psi0(:, psi_index) = temp;
    end
    
elseif strcmp( psi_type, 'rand_pro')
    psi0 = zeros(2^N_atom , psi_N);
    for psi_index = 1: psi_N
        temp = 1;
        for N_index = 1:N_atom
            this_spin = 2* rand(2,1) -1 + i * (2* rand(2,1) -1);
            this_spin = this_spin / norm(this_spin);
            temp = kron(temp, this_spin);
        end
        psi0(:, psi_index) = temp;
    end
else
    error('no this psi0_type')
end
% if strcmp(psi0_type, 'rand')

%psi0_t = psi0';

end

