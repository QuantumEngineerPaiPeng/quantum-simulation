%20190428
% compute observable ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty,
% with disorder
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 4/28

%% runtime parameters
N_list = [4:6];

hx = 1;
Jz = 1;
hz = 0;
repeat = 20;
hz_dis = 2;
hz_dis_list = [repelem(hz_dis, repeat)];

hdtau = 0.1;
htau_list = [hdtau:hdtau:3];

t_list = 20;
t_intervals = {[t_list, t_list], -100}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% begin of for
tic
[ZZ_list, XX_list] = deal( zeros( length(N_list),length(hz_dis_list), length(t_intervals), length(htau_list)) );
for N_index = 1: length(N_list)
    N_atom = N_list(N_index)
    for hz_index = 1: length(hz_dis_list)
        hz_dis = hz_dis_list(N_index)
        [ZZ_list(N_index, hz_index, :,:), XX_list(N_index, hz_index, :,:)] = dipolar_ZZ( ...,
            N_atom, hx, Jz, hz ...
            ,hdtau, htau_list, t_intervals, {{'load_save', load_save}, {'ising', ising}, {'hz_dis', hz_dis}});
    end
end
toc

%% average and variance of different realizations
[XX_mean, XX_sqr_var] = deal( zeros( length(N_list), length(t_intervals), length(htau_list)) );
for N_index = 1: length(N_list)
    N_atom = N_list(N_index);
    XX_mean(N_index, :, :) = mean(XX_list( N_index,:, :, :), 2 );
    XX_sqr_var(N_index, :, :) = sqrt( var(XX_list( N_index,:, :, :), [], 2 ) /repeat );
end

%% derivative, sharpness
%{
figure(49)
clf;

subplot(1,2,1)
normal_factor = 1;% 1./XX_mean(:,1,1) .* XX_mean(1,1,1);
diff_ave = squeeze( diff( XX_mean(:,1,:).* normal_factor, [], 3)  ) ;
diff_err = squeeze( (XX_sqr_var(:,1,4:end) + XX_sqr_var(:,1,1:end-3)).* normal_factor);
ycPlot(htau_list(1:end-3), diff_ave(:,1:end-2) + diff_ave(:,2:end-1) + diff_ave(:,3:end) ...
    , '$h_x\tau$', '$\partial_\tau <X(t)X(0)>_{\beta=0}$', ...,
    ['$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=',mat2str(t_list), ...
    ' ,N_{realize}=',mat2str(repeat),'$'], legends, {{'err', diff_err }});

subplot(1,2,2)
normal_factor = 1./XX_mean(:,2,1) .* XX_mean(1,2,1);
diff_ave = squeeze( diff( XX_mean(:,2,:).* normal_factor, [], 3)  ) ;
diff_err = squeeze( (XX_sqr_var(:,2,4:end) + XX_sqr_var(:,2,1:end-3)).* normal_factor);
ycPlot(htau_list(1:end-3), diff_ave(:,1:end-2) + diff_ave(:,2:end-1) + diff_ave(:,3:end) ...
    , '$h_x\tau$', '$\partial_\tau <X(t)X(0)>_{\beta=0}$', ...,
    ['$J_zt=\infty$, normalized by $XX|_{\tau=0}$'], legends, {{'err', diff_err }});
%}

%% plot average and errorbar
legends = cell(0);
for N_index = 1:length(N_list)
    legends{end+1} = ['$N_{atom}=',mat2str(N_list(N_index)),'$'];
end

figure(48);
clf;

subplot(1,2,1);
normal_factor = 1;
ycPlot(htau_list, squeeze( XX_mean(:,1,:).* normal_factor...%
    ), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['$\frac{h_x}{J_z}=',mat2str(hx), '\frac{h^{dis}_z}{J_z}=',mat2str(hz_dis) ...
    , ' ,J_zt=',mat2str(t_list), ...
    ' ,N_{realize}=',mat2str(repeat),'$'], legends, {{'err', squeeze( XX_sqr_var(:,1,:).* normal_factor)}});

subplot(1,2,2);
normal_factor = 1; % ./XX_mean(:,2,1) .* XX_mean(1,2,1);
ycPlot(htau_list', squeeze( XX_mean(:,2,:) .* normal_factor ...
    ), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['$\frac{h_x}{J_z}=',mat2str(hx) ...
    ,' ,J_zt=\infty$'], legends, {{'err', squeeze( XX_sqr_var(:,2,:).* normal_factor)}});

profile viewer



