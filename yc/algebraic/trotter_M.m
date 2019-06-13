%20190313
% compute observable M of trotterized Ising chain with algebraically decaying
% interacions
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 8;

alpha = 100;
hx = 1/4;
Jz = 1;

Jza = Jz * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [0, 0, 1];
bc = 'o';
sym = [];

dtau = 0.2;
tau_list = [dtau:dtau:4];

psi_N = 4;
psi_type = 'up/down'; % 'up': all up, 'up/down': each spin either up or down, product state,
     % 'rand': total random in 2^N Hilbert space, 'rand_pro': each spin
     % rand direction, different spin no entangle

F_N = 10^3;

filename=sprintf('E:/career/junior2/code/yc/algebraic/alpha%dN%dhx025dtau02.mat',alpha,N_atom);
     
profile on
profile clear

%% test state of IPR
psi0 = initial_psi(psi_type, N_atom, psi_N);
psi0_t = psi0';

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

if exist(filename,'file')
    load(filename);
    loaded = true;
else
    Ut_all = cell(0);
    Ut=OperatorClass(N_atom);
    Ut.matrix = { eye(2^N_atom) };
    U1 = expm( -i* Hz.matrix{1}* dtau );
    U2 = expm( -i* Hx.matrix{1}* dtau );
    
    loaded = false;
end % if exist(filename,'file')

%% test operator 
M0 = 1/N_atom * OperatorClass(N_atom,'z',1/2).matrix{1};

%% begin of for
M_list = zeros(psi_N, length(tau_list));

for tau_index = 1:length(tau_list)

tau = tau_list(tau_index)
if loaded
    Ut = Ut_all{tau_index};
else
    Ut.matrix{1} = U1 * Ut.matrix{1} * U2;
    Ut_all{end+1} = copy( Ut);
end

%% compute M(t -> \infty)
tic;
    psi = Ut.matrix{1}^F_N * psi0;
    temp = 0;
    for t = 1:F_N
        psi = Ut.matrix{1} * psi;
        temp = temp + diag( psi' * M0 * psi);  % W0' = W0
    end % for t = 1:IPR_N

    M_list(:, tau_index) = real( temp) / F_N;

toc;
    
end % for len_index = 1:length(N_atom_list)

%% after for, save
if ~loaded
    save(filename, 'Ut_all');
end

%% plot log_IPR vs N_atom
legends = cell(0);
for psi_index = 1:psi_N
    legends{end+1} = ['$N_{up}=',mat2str(N_atom +1 -psi_index),'$'];
end

ycPlot(14, tau_list, M_list, '$J_z\tau$', '$M$', ...,
    ['$\alpha=',mat2str(alpha),' ,h_x/J_z=',mat2str(hx),' ,N_{atom}=',mat2str(N_atom), ...
    ' ,N_{period}=',mat2str(F_N),'$'], legends);

profile viewer