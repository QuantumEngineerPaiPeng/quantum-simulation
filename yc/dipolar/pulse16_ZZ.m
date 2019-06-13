%2019 \pi
% compute observable ZZ(t->\infty) and XX in experiment pulses with
% different \tau, dipolar interaction at T = \infty
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U(t) Z_0 ) / 2^N, Z_0 = S^z
% H_{\tau=0}: dipolar_y 
% = Jz \sum_{i<j} (S^y_i S^y_j -1/2 S^x_i S^x_j -1/2 S^z_i S^z_j) / |i-j|^alpha
%
% \tau \neq 0 pulse: t1 x t2 y 2t1 y t2 x, 2t1 x t2 y 2t1 y t2 x, 2t1 -x t2
% -y 2t1 -y t2 -x, 2t1 -x t2 -y 2t1 -y t2 -x t1,
% x:=e^{-ipi X}, t1 = tau(1-u), t2 = tau(1+2u), 4t1+2t2=6tau
% Written using OperatorClass
% Created by Chao Yin 6/9


%% runtime parameters
N_all = [6:2:12];

Jz = 1;
u = 0.8;
u_list = [2*u, 4*u, 3*u];
bc = 'p';
sym = ['kp'];
hz = 0;

dtau = 0.05;
tau_list = [dtau:dtau:2];

t_ave = 10^3;
t_interval = [t_ave, t_ave+100];

profile on
profile clear

ZZ_list = zeros(length(N_all), 2, length(tau_list)+1);
for N_index= 1: length(N_all)

N_atom = N_all(N_index)

%% aim hamiltonian and its longtime ZZ
alph = 3;
cplist = Jz * [1:N_atom-1].^(-alph) ;
dipolar_y = [-1/2 1 -1/2];
H_aim = u* OperatorClass(N_atom,dipolar_y,1/4,bc,cplist) + hz * OperatorClass(N_atom,'z',1/2);
H_aim.symmetrize(sym);
Y0 = OperatorClass(N_atom,'y',1/2);
Y0.symmetrize(sym);
d0 = H_aim* Y0 - Y0* H_aim;
U_aim = H2U(H_aim, 1);

%% test operator 
Z0 = {Y0, 1/u* H_aim};
for oper_index = 3: 1: length(Z0)
    Z0{oper_index}.symmetrize(sym);
end

ZZ_list(N_index, :, 1) = ZZ_longtime( U_aim, t_interval, Z0)*4/N_atom;

%% begin of for
for tau_index = 1:length(tau_list)

tau = tau_list(tau_index);
tau_whole = 24* tau;

Ut = U_16pulse(N_atom, bc, sym, Jz, u, tau, hz );

ZZ_list(N_index, :,tau_index+1) = ZZ_longtime(Ut, t_interval, Z0)*4/N_atom;

end % for len_index = 1:length(N_atom_list)
end % for N_index

%% plot log_IPR vs N_atom
figure(41)
clf;
legends = {};
for N_index = 1: length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str(N_all(N_index)),'$'];
end


subplot(1,2,1)
ycPlot([0, tau_list], squeeze( ZZ_list(:,1,:)), '$J_z\tau$', '$<Y(\infty)Y(0)>_{\beta=0}$', ...,
    ['Ken16 $Dip_Y, h_z=',mat2str(hz),' ,u=',mat2str(u), ' ,n_{period}=1e3:1e3+100', ...
    '$'], legends, {{'x_par', pi/24/hz}});
subplot(1,2,2)

ycPlot([0, tau_list], squeeze( ZZ_list(:,2,:)), '$J_z\tau$', '$<Dip_Y(\infty)Dip_Y(0)>_{\beta=0}$', ...,
    [], legends, {{'x_par', pi/24/hz}});


profile viewer

%{
validity plot
ycPlot(40, tau_list, delta_max, '$J_z\tau$', '$|U_{real}-U_{aim}|$' ...,
    , ['$N_{atom}=', mat2str(N_atom),'$']);
%}




