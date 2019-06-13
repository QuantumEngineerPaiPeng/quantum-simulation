%2019 \pi
% confirm the 16 pulse is valid
% H_{\tau=0}: dipolar_y 
% = Jz \sum_{i<j} (S^y_i S^y_j -1/2 S^x_i S^x_j -1/2 S^z_i S^z_j) / |i-j|^alpha
%
% \tau \neq 0 pulse: t1 x t2 y 2t1 y t2 x, 2t1 x t2 y 2t1 y t2 x, 2t1 -x t2
% -y 2t1 -y t2 -x, 2t1 -x t2 -y 2t1 -y t2 -x t1,
% x:=e^{-ipi X}, t1 = tau(1-u), t2 = tau(1+2u), 4t1+2t2=6tau
% Written using OperatorClass
% Created by Chao Yin


%% runtime parameters
N_atom = 10;

Jz = 1;
u = 0.2;
bc = 'o';
sym = [];
hz = 0;

tau_list = 10.^[-5:1];
t_test = 10^3;
profile on
profile clear

%% begin of for
dist = tau_list;
for tau_index = 1: length(tau_list)

tau = tau_list(tau_index);
tau_whole = 24*tau;
period_N = floor( t_test/ tau_whole );

%% aim hamiltonian and its Floquet
alph = 3;
cplist = Jz * [1:N_atom-1].^(-alph) ;
dipolar_y = [-1/2 1 -1/2];
H_aim = u* OperatorClass(N_atom,dipolar_y,1/4,bc,cplist) + hz * OperatorClass(N_atom,'z',1/2);
U_aim = H2U(H_aim, tau_whole * period_N);

%% pulse Floquet
U_real = U_16pulse(N_atom, bc, sym, Jz, u, tau, hz);

%% distance
dU = U_real^period_N - U_aim;
dist(tau_index) = max(max( abs(dU.matrix{1}) ));

end % for tau_index = 1: length(tau_list)

%% plot

ycPlot(40, tau_list, dist, '$J_z\tau$', '$|U_{real}-U_{aim}|$' ...,
    , ['$N_{atom}=', mat2str(N_atom),' ,T_{test}=', mat2str(t_test),'$']);

profile viewer
%}


