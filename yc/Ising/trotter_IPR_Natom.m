%20190302
% compute IPR of trotterized traverse Ising chain
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom_list = [13];

model = [0, 0, 1];
bc = 'o';
cplist = [1];
sym = [];

g = 2;
h = 2;

IPR_list = N_atom_list;

%% begin of for
for len_index = 1:length(N_atom_list)

N_atom = N_atom_list(len_index);

%% test state
psi0 = zeros(2^N_atom ,1);
psi0(1) = 1;
psi0_t = psi0';

%% construct Hamiltonian before trotter
H = OperatorClass(N_atom,model,1/4,bc,cplist) + h* OperatorClass(N_atom,'z',1/2) ...
    + g* OperatorClass(N_atom,'x',1/2);

%% diagonalize
tic;
H.diagonalize();
toc;

%% IPR of psi0
IPR_list(len_index) = sum( abs( psi0_t * H.eigsys{1}.V).^4 );

end % for len_index = 1:length(N_atom_list)

%% after for

IPR_list
log_IPR = -log2(IPR_list) ./ ( N_atom_list - 1 )

%% plot log_IPR vs N_atom
ycPlot(1, N_atom_list, log_IPR, '$N_{atom}$', '$\lambda_{IPR} / \lambda_{D}$', ...,
    ['$g=',mat2str(g),' ,h=',mat2str(h),' ,\tau=0$']);








