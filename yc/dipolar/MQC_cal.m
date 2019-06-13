%20190323
% compute MQC of algebraic Ising with hx and hz
% interacions at T = \infty
% S_mt = S_m(t) = 2^(-L) tr( O0 U(m,t)' O0 U(m,t) )
% U(m, t) = U_t^(-t) exp(im pi PZ/L) U_t^t
%
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/23

%% runtime parameters
N_atom = 3;

hx = 1;
Jz = 1;
hz = 0;

hdtau = pi/2;
htau_list = [hdtau];

t_start = 8;
t_end = t_start;
t_intervals = {[t_start, t_end]}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% use function MQC
[I_qtz, I_qtx] = MQC( N_atom, hx, Jz, hz ...
    ,hdtau./ hx, htau_list./ hx, t_intervals, load_save, ising);

%% fix tau
figure(53);
clf;
curve_n = 7;
tau_index = 157;
legends = cell(0);
for q_index = 0: N_atom
    legends{end+1} = ['$q=',mat2str(q_index),'$'];
end
x_data = linspace(t_intervals{1}(1), t_intervals{1}(2) , size( I_qtz{1,tau_index},2) );

subplot(1,2,1);
ycPlot(x_data, I_qtz{1,tau_index}(1: curve_n , :), '$J_z t$', '$I_q(Z;t)$', ...,
    ['dipolar z, $\frac{h_z}{J_z}=',mat2str(hz),' ,\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,h_x \tau=', mat2str(htau_list(tau_index)) ,'$'], legends);
ylim([0 1])
subplot(1,2,2);
ycPlot(x_data, I_qtx{1,tau_index}(1: curve_n , :), '$J_z t$', '$I_q(X;t)$', ...,
    ['dipolar z,$\frac{h_z}{J_z}=',mat2str(hz),' ,\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,h_x \tau=', mat2str(htau_list(tau_index)) ,'$'], legends);

%% x_data: tau, fix t=last
figure(54);
clf;
y_data = zeros( curve_n, length(htau_list), 2);
for tau_index = 1: length(htau_list)
    y_data(:, tau_index, 1) = mean( I_qtz{1,tau_index}(1:curve_n , :), 2 );
    y_data(:, tau_index, 2) = mean( I_qtx{1,tau_index}(1:curve_n , :), 2 );
end

subplot(1,2,1);
ycPlot(htau_list, y_data(:,:,1), '$h_x\tau$', '$I_q(Z;t)$', ...,
    ['dipolar z,$\frac{h_z}{J_z}=',mat2str(hz),' ,\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=', mat2str(t_end) ,'$'], legends);
ylim([0 1])
subplot(1,2,2);
ycPlot(htau_list, y_data(:,:,2), '$h_x\tau$', '$I_q(X;t)$', ...,
    ['dipolar z,$\frac{h_z}{J_z}=',mat2str(hz),' ,\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=', mat2str(t_end) ,'$'], legends);


profile viewer







