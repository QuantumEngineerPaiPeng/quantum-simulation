%20190319
% compute observable ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty, change hx while hz=0
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/28

%% runtime parameters
N_atom = 8;

hx_list = [0.2, 1, 1.5];
Jz = 1;
hz = 0;

hdtau = 0.1;
htau_list = [hdtau:hdtau:3];

t_list = 50;
t_intervals = {-100}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% begin of for
[ZZ_list, XX_list] = deal( zeros(length(hx_list), length(t_intervals), length(htau_list)) );
for hx_index = 1: length(hx_list)
    hx = hx_list(hx_index);
    [ZZ_list(hx_index, :,:), XX_list(hx_index, :,:)] = dipolar_ZZ( ...,
        N_atom, hx, Jz, hz ...
        ,hdtau./hx, htau_list./hx, t_intervals, {{'load_save',load_save}} );
end

%% exact result of ising
cri_htau = pi ./ (1+ 1./ (2* hx_list) );

%% plot ZZ, XX vs N_atom

legends = cell(0);
for hx_index = 1:length(hx_list)
    legends{end+1} = ['$\frac{h_x}{J_z}=',mat2str(hx_list(hx_index)),'$'];
end

figure(47);
clf;

subplot(1,2,1);
ycPlot(htau_list, squeeze( XX_list(:,1,:)), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=',mat2str(t_list),'$'], legends);
hold on
subplot(1,2,2);
ycPlot(htau_list, squeeze( ZZ_list(:,1,:)), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);

%% ising case
figure(48);
clf;

%subplot(1,2,1);
extra = cell(0);
for hx_index = 1: length(hx_list)
    extra{end+1} = {'x_par', cri_htau(hx_index)};
end
ycPlot(htau_list, squeeze(  XX_list(:,2,:)), '$h_x\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends, extra);


subplot(1,2,2);
ycPlot(htau_list, squeeze( ZZ_list(:,2,:)), '$h_x\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['ising,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$'], legends);
%}


%% tau \sim 0, == idea traverse ising model
figure(49);
clf;

%subplot(1,2,1);
ycPlot(hx_list, squeeze(  XX_list(:,1,:)), '$h_x/J_z$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['ising, $J_z\tau=', mat2str(hdtau),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,J_zt=\infty$']);


profile viewer



