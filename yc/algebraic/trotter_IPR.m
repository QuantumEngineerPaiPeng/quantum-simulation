%20190311
% compute IPR of trotterized Ising chain with algebraically decaying
% interacions
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 10;

alpha = 0;
hx = 1;
Jz = sqrt(17);

Jza = Jz * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [0, 0, 1];
bc = 'o';
sym = [];

IPR_ED = true;
IPR_N = 10^3;

dtau = 0.2;
tau_list = [dtau:dtau:4];

psi_N = 1;
psi_type = 'up'; % 'up': all up, 'up/down': each spin either up or down, product state,
     % 'rand': total random in 2^N Hilbert space, 'rand_pro': each spin
     % rand direction, different spin no entangle

profile on
profile clear

%% test state of IPR


psi0_t = psi0';

%% construct Hamiltonian before trotter
Hz=OperatorClass(N_atom,model,1/4,bc,cplist);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

Ut=OperatorClass(N_atom);
Ut.matrix = { eye(2^N_atom) };

%% U1, U2 for dtau
U1 = expm( -i* Hz.matrix{1}* dtau );
U2 = expm( -i* Hx.matrix{1}* dtau );

%% begin of for
IPR_list = zeros(psi_N, length(tau_list));

for tau_index = 1:length(tau_list)

tau = tau_list(tau_index)

tic;
Ut.matrix{1} = U1 * Ut.matrix{1} * U2;
toc;

%% IPR of psi0
tic;
if IPR_ED == true
    Ut.diagonalize();
    %temp = sum( abs( psi0_t * Ut.eigsys{1}.V).^4 , 2);
    IPR_list(:, tau_index) = sum( abs( psi0_t * Ut.eigsys{1}.V).^4 , 2);
else
    psi = Ut.matrix{1}^IPR_N * psi0;
    temp = 0;
    for t = 1:IPR_N
        psi = Ut.matrix{1} * psi;
        temp = temp + abs( psi0_t * psi )^2;
    end % for t = 1:IPR_N

    IPR_list(tau_index) = temp / IPR_N;
end
toc;
    
end % for len_index = 1:length(N_atom_list)

%% after for
log_IPR = -log2(IPR_list) ./ ( N_atom - 3 );
log_IPR_all(end+1, :) = log_IPR;
J_list(end+1) = Jz;

%% plot log_IPR vs N_atom
%{
legends = cell(0);
for psi_index = 2:length(N_all)
    legends{end+1} = ['$N=',mat2str(N_all(psi_index)),'$'];
end
%}

ycPlot(11, Jz* tau_list, log_IPR, '$J_z\tau$', '$\lambda_{IPR} / \lambda_{D}$', ...,
    ['$\alpha=', mat2str(alpha),' ,h_x=',mat2str(hx),' ,J_z=',mat2str(Jz), ' ,N_{atom}=',mat2str(N_atom),'$'] ...
    );

profile viewer
