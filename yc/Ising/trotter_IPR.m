%20190302
% compute IPR of trotterized traverse Ising chain
% H = HZ + Hx
% HZ = \sum S^z_j S^z_{j+1} + h \sum S^z_j  (J=1)
% Hx = g \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 8;

model = [0, 0, 1];
bc = 'o';
cplist = [1];
sym = [];

g = 1;
h = 1;

IPR_ED = true;
IPR_N = 30;

dtau = 0.2;
tau_list = [dtau:dtau:2];

psi_N = 2;
psi_type = 'up/down'; % 'up': all up, 'up/down': each spin either up or down, product state,
     % 'rand': total random in 2^N Hilbert space, 'rand_pro': each spin
     % rand direction, different spin no entangle

profile on
profile clear

%% test state of IPR
psi0 = initial_psi(psi_type, N_atom, psi_N);
psi0_t = psi0';

%% construct Hamiltonian before trotter
H_int=OperatorClass(N_atom,model,1/4,bc,cplist);  % H_int = \sum \sigma^z_j \sigma^z_{j+1}
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1/2);  % Hz = \sum \sigma^z_j
Hz.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

HZ = H_int + h*Hz;

H = HZ + g*Hx;

Ut=OperatorClass(N_atom);
Ut.matrix = { eye(2^N_atom) };

%Ut_all = cell(10);

%% U1, U2 for dtau
U1 = expm( -i* HZ.matrix{1}* dtau );
U2 = expm( -i*g * Hx.matrix{1}* dtau );

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
%log_IPR_all(end+1, :) = log_IPR;
%g_list(end+1) = g;

%% plot log_IPR vs N_atom
legends = cell(0);
for psi_index = 1:psi_N
    legends{end+1} = ['$N_{up}=',mat2str(0),'$'];
end

ycPlot(4, tau_list, log_IPR, '$J\tau$', '$\lambda_{IPR} / \lambda_{D}$', ...,
    ['$g=',mat2str(g),' ,h=',mat2str(h),' ,N_{atom}=',mat2str(N_atom),'$'], ...
    legends);
%xticks([0:2*pi:10*pi])
%xticklabels({'0','2pi','4\pi','6\pi','8\pi'})

profile viewer
