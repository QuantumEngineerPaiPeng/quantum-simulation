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
N_atom = 8;

hx = 1;
Jz = 1;
hz = 0;
repeat = 20;
hz_dis_list = [0, repelem(0.5, repeat), repelem(1, repeat), repelem(2, repeat)];

hdtau = 0.05;
htau_list = [hdtau:hdtau:3];

t_list = 20;
t_intervals = {[t_list, t_list], -100}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% begin of for
tic
[ZZ_list, XX_list] = deal( zeros(length(hz_dis_list), length(t_intervals), length(htau_list)) );
for dis_index = 1: length(hz_dis_list)
    hz_dis = hz_dis_list(dis_index)
    [XX_list(dis_index, :,:), ZZ_list(dis_index, :,:)] = dipolar_ZZ( ...,
        N_atom, hx, Jz, hz ...
        ,hdtau, htau_list, t_intervals, {{'load_save', load_save}, {'ising', ising}, {'hz_dis', hz_dis}});
end
toc

%% plot ZZ, XX vs N_atom

legends = cell(0);
for dis_index = 1:length(hz_dis_list)
    legends{end+1} = ['$\frac{h^{dis}_z}{J_z}=',mat2str(hz_dis_list(dis_index)),'$'];
end

figure(47);
clf;

y_label = '$<H_{dipz}(t)H_{dipz}(0)>_{\beta=0}$';
subplot(1,2,1);
ycPlot(htau_list, squeeze( XX_list(:,1,:)), '$h_x\tau$', y_label, ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=',mat2str(t_list),'$'], legends);
hold on
subplot(1,2,2);
ycPlot(htau_list, squeeze( XX_list(:,2,:)), '$h_x\tau$', y_label,...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);

%% average and variance of different realizations
XX_mean = XX_list(1, :, :);
XX_sqr_var = zeros(1, 2, length(htau_list));
for dis_index = 1:3
    XX_mean(end+1, :, :) = mean(XX_list( (dis_index-1)*repeat+2 : dis_index*repeat+1, :, :), 1 );
    XX_sqr_var(end+1, :, :) = sqrt( var(XX_list( (dis_index-1)*repeat+2 : dis_index*repeat+1, :, :), [], 1 ) /repeat );
end % for dis_index = 1:3

%% derivative, sharpness
figure(49)
clf;

subplot(1,2,1)
normal_factor = 1./XX_mean(:,1,1) .* XX_mean(1,1,1);
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

%% 
legends = cell(0);
hz_mean_list = [0 0.5 1 2];
for dis_index = 1:length(hz_mean_list)
    legends{end+1} = ['$\frac{h^{dis}_z}{J_z}=',mat2str(hz_mean_list(dis_index)),'$'];
end

figure(48);
clf;

subplot(1,2,1);
normal_factor = 1;% 1./XX_mean(:,1,1) .* XX_mean(1,1,1);
ycPlot(htau_list, squeeze( XX_mean(:,1,:).* normal_factor...%
    ), '$h_x\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=',mat2str(t_list), ...
    ' ,N_{realize}=',mat2str(repeat),'$'], legends, {{'err', squeeze( XX_sqr_var(:,1,:).* normal_factor)}});

subplot(1,2,2);
normal_factor = 1;% 1./XX_mean(:,2,1) .* XX_mean(1,2,1);
ycPlot(htau_list, squeeze( XX_mean(:,2,:) .* normal_factor ...
    ), '$h_x\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['$\frac{h_x}{J_z}=',mat2str(hx) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends, {{'err', squeeze( XX_sqr_var(:,2,:).* normal_factor)}});

profile viewer

%%
%{
figure(21)
ycPlot(htau_list.^2, squeeze( XX_list(:,2,:)), '$(h_x\tau)^2$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);
%}



