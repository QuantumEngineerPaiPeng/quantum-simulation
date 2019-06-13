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
N_atom_list = [4: 12];

alpha = 0;
hx = 1;
Jz = 4;

cplist = cell(0); % Jza = cplist{i}[1]
for N_index = 1: length(N_atom_list)
    N_atom = N_atom_list(N_index);
    Jza = Jz * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^(-alpha)  ) );
    cplist{end+1} = Jza * [1:N_atom-1].^(-alpha) ;
end

model = [0, 0, 1];
bc = 'o';
sym = [];

IPR_list = N_atom_list;

%% begin of for
for len_index = 1:length(N_atom_list)

N_atom = N_atom_list(len_index);

%% test state
psi0 = zeros(2^N_atom ,1);
psi0(1) = 1;
psi0_t = psi0';

%% construct Hamiltonian before trotter
H = OperatorClass(N_atom,model,1/4,bc,cplist{len_index}) + hx * OperatorClass(N_atom,'x',1/2);

%% diagonalize
tic;
H.diagonalize();
toc;

%% IPR of psi0
IPR_list(len_index) = sum( abs( psi0_t * H.eigsys{1}.V).^4 );

end % for len_index = 1:length(N_atom_list)

%% after for

log_IPR = -log2(IPR_list) ./ ( N_atom_list - 3 ); % \lambda_D = log(2)(1-3/N)

%% plot log_IPR vs N_atom
ycPlot(1, N_atom_list, log_IPR, '$N_{atom}$', '$\lambda_{IPR} / \lambda_{D}$', ...,
    ['$\alpha=',mat2str(alpha),' ,h_x=',mat2str(hx), ' ,J_z=',mat2str(Jz), ' ,\tau=0$']);








