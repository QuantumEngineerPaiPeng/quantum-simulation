%20190303
% compute OTO correlator F of trotterized traverse Ising chain
% H = HZ + Hx
% HZ = \sum S^z_j S^z_{j+1} + h \sum S^z_j
% Hx = g \sum S^x_j
% F(t) = < V^d(t) W^d V(t) W >, V = W = N^{-1} \sum S^z_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 8;

model = [0, 0, 1];
bc = 'o';
cplist = [1];
sym = [];

g = 2;
h = 2;

dtau = 0.1;
tau_list = [dtau:dtau:1];

F_N = 10^4;

%% test state of OTOC
psi0 = zeros(2^N_atom ,1);
psi0(1) = 1;
psi0_t = psi0';

%% construct Hamiltonian before trotter
H_int=OperatorClass(N_atom,model,1/4,bc,cplist);
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

HZ = H_int + h*Hz;

Ut=OperatorClass(N_atom);
Ut.matrix = { eye(2^N_atom) };

%% test operator of OTOC
W0 = 1/N_atom * Hz.matrix{1};

%% U1, U2 for dtau
U1 = expm( -i* HZ.matrix{1}* dtau );
U2 = expm( -i*g * Hx.matrix{1}* dtau );

%% begin of for
F_list = tau_list;

for tau_index = 1:length(tau_list)

tau = tau_list(tau_index)
tic;
Ut.matrix{1} = U1 * Ut.matrix{1} * U2;
toc;

%% compute F(t -> \infty)
tic;
    Ut_FN = Ut.matrix{1}^F_N;
    Vt = Ut_FN' * W0 * Ut_FN;
    temp = 0;
    for t = 1:F_N
        Vt = Ut.matrix{1}' * Vt * Ut.matrix{1};
        temp = temp + psi0_t * Vt' * W0 * Vt * W0 * psi0;  % W0' = W0
    end % for t = 1:IPR_N

    F_list(tau_index) = real( temp) / F_N;

toc;
    
end % for len_index = 1:length(N_atom_list)

%% after for

%% plot log_IPR vs N_atom
ycPlot(3, tau_list, 8*F_list, '$\tau$', '$F / F_0$', ...,
    ['$g=',mat2str(g),' ,h=',mat2str(h),' ,N_{atom}=',mat2str(N_atom),' ,N_{period}=',mat2str(F_N),'$']);


