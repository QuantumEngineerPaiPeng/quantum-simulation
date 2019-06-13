%20190515
% compute observable XX(t) of ising with hx 
% at T = \infty, change interaction decay rate alpha, hz=0
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/15

%% runtime parameters
N_atom = 12;

hx = 0.5;
Jz = 1;
hz = 0;

alpha_list = [0 100];

hdtau = 0.05;
htau_list = [hdtau:hdtau:3];

t_list = 50;
t_intervals = {-100}; % negative number triggers ED for infinite long time

load_save = false;
ising = true;

profile on
profile clear

%% begin of for
[ZZ_list, XX_list] = deal( zeros(length(alpha_list), length(t_intervals), length(htau_list)) );
for a_index = 1: length(alpha_list)
    alpha = alpha_list(a_index);
    [ZZ_list(a_index, :,:), XX_list(a_index, :,:)] = dipolar_ZZ( ...,
        N_atom, hx, Jz, hz ...
        ,hdtau /hx, htau_list/ hx, t_intervals, {{'ising',ising}, {'alpha', alpha}});
end

%% exact result of ising
cri_htau = pi ./ (1+ 1./ (2* hx) );

%% ising case
legends = cell(0);
for a_index = 1:length(alpha_list)
    legends{end+1} = ['$\alpha=',mat2str(alpha_list(a_index), 3),'$'];
end

figure(60);
clf;

%subplot(1,2,1);
extra = { {'x_par', cri_htau} };

ycPlot(htau_list, squeeze(  XX_list(:,1,:)), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends, extra);

%{
subplot(1,2,2);
ycPlot(htau_list, squeeze( ZZ_list(:,1,:)), '$h_x\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);
%}

profile viewer



