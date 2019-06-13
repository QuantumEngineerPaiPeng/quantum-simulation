%20190303
% compute observable M = W of trotterized traverse Ising chain
% H = HZ + Hx
% HZ = \sum S^z_j S^z_{j+1} + h \sum S^z_j
% Hx = g \sum S^x_j
%  W = N^{-1} \sum S^z_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 8;

model = [0, 0, 1];
bc = 'o';
cplist = [1];
sym = [];

g = 0;
h = 1;

dtau = 0.2;
tau_list = [dtau:dtau:4];

psi_N = 4;
psi_type = 'up/down'; % 'up': all up, 'up/down': each spin either up or down, product state,
     % 'rand': total random in 2^N Hilbert space, 'rand_pro': each spin
     % rand direction, different spin no entangle

F_N = 10^3;

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

psi0_t = psi0';

%% construct Hamiltonian before trotter
H_int=OperatorClass(N_atom,model,1/4,bc,cplist);
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

HZ = H_int + h*Hz;

H = HZ + g*Hx;

Ut=OperatorClass(N_atom);
Ut.matrix = { eye(2^N_atom) };

%% test operator of OTOC
W0 = 1/N_atom * Hz.matrix{1};

%% U1, U2 for dtau
U1 = expm( -i* HZ.matrix{1}* dtau );
U2 = expm( -i*g * Hx.matrix{1}* dtau );

%% begin of for
M_list = zeros(psi_N, length(tau_list));

for tau_index = 1:length(tau_list)

tau = tau_list(tau_index)
tic;
Ut.matrix{1} = U1 * Ut.matrix{1} * U2;
toc;

%% compute M(t -> \infty)
tic;
    psi = Ut.matrix{1}^F_N * psi0;
    temp = 0;
    for t = 1:F_N
        psi = Ut.matrix{1} * psi;
        temp = temp + diag( psi' * W0 * psi );  % W0' = W0
    end % for t = 1:IPR_N

    M_list(:, tau_index) = real( temp) / F_N;

toc;
    
end % for len_index = 1:length(N_atom_list)

%% after for

%% plot log_IPR vs N_atom
legends = cell(0);
for psi_index = 1:psi_N
    legends{end+1} = ['$N_{up}=',mat2str(N_atom +1 -psi_index),'$'];
end

ycPlot(4, tau_list, M_list, '$J\tau$', '$M$', ...,
    ['$g=',mat2str(g),' ,h=',mat2str(h),' ,N_{atom}=',mat2str(N_atom),' ,N_{period}=',mat2str(F_N),'$'], ...
    legends);

