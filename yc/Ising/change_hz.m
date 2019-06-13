%20190329
% compute observable ZZ(t) and XX(t) of ising with hx and hz
% interacions at T = \infty, change hz while hx=1, from hz=0 (solvable) to
% hz=1 (heyl's paper)
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/29

%% runtime parameters
N_atom = 6;

hx = 1;
Jz = 1;
hz_list = [0:0.2:1];

hdtau = 0.05;
htau_list = [hdtau:hdtau:3];

t_list = 50;
t_intervals = {-100}; % negative number triggers ED for infinite long time

load_save = false;
ising = true;

profile on
profile clear

%% begin of for
[ZZ_list, XX_list] = deal( zeros(length(hz_list), length(t_intervals), length(htau_list)) );
for hz_index = 1: length(hz_list)
    hz = hz_list(hz_index);
    [ZZ_list(hz_index, :,:), XX_list(hz_index, :,:)] = dipolar_ZZ( ...,
        N_atom, hx, Jz, hz ...
        ,hdtau /hx, htau_list/ hx, t_intervals, {{'ising',true}});
end

%% exact result of ising
cri_htau = pi ./ (1+ 1./ (2* hx) );

%% ising case
legends = cell(0);
for hz_index = 1:length(hz_list)
    legends{end+1} = ['$\frac{h_z}{J_z}=',mat2str(hz_list(hz_index)),'$'];
end

figure(48);
clf;

subplot(1,2,1);
extra = { {'x_par', cri_htau} };

ycPlot(htau_list, squeeze(  XX_list(:,1,:)), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends, extra);


subplot(1,2,2);
ycPlot(htau_list, squeeze( ZZ_list(:,1,:)), '$h_x\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);
%}

profile viewer



