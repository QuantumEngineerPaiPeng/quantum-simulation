%20190509
% compute observable <dipo(t)dipo(0)> of dipolar_z with hx and hz
% interacions at T = \infty,
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/9

%% runtime parameters
N_atom_list = [6:2:12];

hx = 0.5;
Jz = 1;
hz = 0;
hz_dis_list = [0];

hdtau = 0.05;
htau_list = [hdtau:hdtau:4];

t_list = 20;
t_intervals = {[t_list, t_list], -100}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% begin of for
tic
[ZZ_list, XX_list] = deal( zeros(length(N_atom_list), length(t_intervals), length(htau_list)) );
for hx_index = 1: length(N_atom_list)
    N_atom = N_atom_list(hx_index)
    [ZZ_list(hx_index, :,:), XX_list(hx_index, :,:)] = dipolar_ZZ( ...,
        N_atom, hx, Jz, hz ...
        ,hdtau/hx, htau_list/hx, t_intervals, {{'load_save', load_save}, {'ising', ising}});
end
toc

%% plot ZZ, XX vs N_atom

legends = cell(0);
for dis_index = 1:length(N_atom_list)
    legends{end+1} = ['$N_{atom}=',mat2str(N_atom_list(dis_index)),'$'];
end

figure(47);
clf;

y_label = '$<H_{dipz}(t)H_{dipz}(0)>_{\beta=0}$';
subplot(1,2,1);
ycPlot(htau_list, squeeze( ZZ_list(:,1,:)), '$h_x\tau$', y_label, ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,\frac{h_x}{J_z}=',mat2str(hx),' ,J_zt=',mat2str(t_list),'$'], legends);
hold on
subplot(1,2,2);
ycPlot(htau_list, squeeze( ZZ_list(:,2,:)), '$h_x\tau$', y_label,...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,\frac{h_x}{J_z}=',mat2str(hx),' ,J_zt=\infty$'], legends);

profile viewer

%%
%{
figure(21)
ycPlot(htau_list.^2, squeeze( XX_list(:,2,:)), '$(h_x\tau)^2$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);
%}



