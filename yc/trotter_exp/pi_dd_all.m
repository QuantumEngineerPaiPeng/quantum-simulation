%20190507
% try to produce the dipolar longtime XX plot from experiment data with DD
% Created by Chao Yin
% last update: 5/01

%% runtime parameters

exp_num_list = [ [20:24, 26], [58: 69] ,[72:83]];
%{
if exp_num_list(1) == 58
    d4 = [2,4,8,16,24,32,40,48,56,64,72,80];
    p2 = 2.95;
elseif exp_num_list(1) == 72
    d4 = [2,4,8,16,24,32,40,48,56,64,72,80];
    p2 = 2.06;
elseif exp_num_list(1) == 20
    d4 = 2.^[1:6];
    p2 = 2.95;
end
%}
d4 = { 2.^[1:6], [2,4,8,16,24,32,40,48,56,64,72,80], [2,4,8,16,24,32,40,48,56,64,72,80] };
p2 = {2.95, 2.95, 2.06};
J = 33e-3;
    
exp_name = 'expt3';

load_data = true;

%% load data
if load_data
    YY_raw = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY_raw{end+1} = B;
    end
end

YY = {[], [], []};
for exp_index = 1: 6
    YY{1}(:, end+1) = imag( YY_raw{exp_index} );
end
for exp_index = 7: 6+12
    YY{2}(:, end+1) = real( YY_raw{exp_index} );
end
for exp_index = 19: 6+12+12
    YY{3}(:, end+1) = real( YY_raw{exp_index} );
end
n_N = size(YY{1},1)

%% fixt
tau_N = length(exp_num_list);
[tau_list, fixt_n_list, fixt_YY, t_list ] = deal( cell(1, 3) );

for run_index = 1: 3
    tau_N = size(YY{run_index}, 2);
    tau_list{run_index} = (d4{run_index}+ p2{run_index});
    t_list{run_index} = [1: n_N]' * tau_list{run_index};
    Jtau_list = J* tau_list{run_index};

    fixt_n_list{run_index} = round( tau_list{run_index}(1)*n_N ./ tau_list{run_index} );
    fixt_YY{run_index} = [];
    for tau_index = 1:tau_N
        fixt_YY{run_index}(end+1) = YY{run_index}(fixt_n_list{run_index}(tau_index), tau_index);
    end
end

%% left: plot all data. right: fix t
figure(13)
clf;
subplot(1,2,1)
legends = {'fix $t$, $\pi_X$', 'fix $t$, $\pi_Y, \tau_1 = 2.95\mu s$', 'fix $t$, $\pi_Y, \tau_1 = 2.06\mu s$'};
for tau_index = 1:tau_N
    %legends{end+1} = ['$\tau=',mat2str( tau_list(tau_index) ),'\mu s$'];
end

hold on
plot_tau = 4;
for run_index = 1:3
    scatter(fixt_n_list{run_index}(1:plot_tau).* tau_list{run_index}(1:plot_tau), fixt_YY{run_index}(1:plot_tau), 'LineWidth', 2);
end

colors = {};
for tau_index = 1: 4
    line = ycPlot(t_list{1}(:, tau_index), YY{1}(:, tau_index), '$t/2/\mu s$', ...
        '$<Y(n\tau)Y(0)>_{\beta=0}$', ['$\pi_Y$-pulse DD'], legends);
    colors{end+1} = get(line, 'color');
    ycPlot(t_list{2}(:, tau_index), YY{2}(:, tau_index), '$t/2/\mu s$', '$<Y(n\tau)Y(0)>_{\beta=0}$'...
        , ['$\pi_Y$-pulse DD'], legends, { {'line_type', '--'},{'color', colors{end}} });
    ycPlot(t_list{3}(:, tau_index), YY{3}(:, tau_index), '$t/2(\mu s)$', '$<Y(n\tau)Y(0)>_{\beta=0}$'...
        , ['$\pi$-pulse DD'], legends, { {'line_type', ':'},{'color', colors{end}} });
end
xlim([0 500])

subplot(1,2,2)
hold on

legends = {'$\pi_X, t=128\times 4.95 \mu s$', '$\pi_Y, \tau_1 = 2.95\mu s, t=128\times 4.95\mu s$'...
    , '$\pi_Y, \tau_1 = 2.06\mu s, t=128\times 4.06\mu s$'};
for run_index = 1:3
    ycPlot(J*tau_list{run_index} * 180/ pi, fixt_YY{run_index}, '$J\tau(^\circ)$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
        ['$\pi$-pulse DD'], legends);
end

%% plot
figure(10)
clf;
subplot(1,2,1)
legends = {'fix $t$'};
for tau_index = 1:tau_N
    legends{end+1} = ['$\tau=',mat2str( tau_list(tau_index) ),'\mu s$'];
end

hold on
scatter(fixt_n_list.* tau_list, fixt_YY, 'LineWidth', 2);
for tau_index = 1: tau_N 
    ycPlot([1: length(YY{tau_index})]* tau_list(tau_index), YY{tau_index}, '$t/2/\mu s$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
        ['$\pi_Y$-pulse DD'], legends);
end
xlim([0 500])

subplot(1,2,2)
hold on

legends =  {'experiment','numeric, $h\tau$'};
ycPlot(Jtau_list * 180/ pi, fixt_YY, '$J\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$t=128\times', mat2str(tau_list(1)),'\mu s$'], legends, {{'Marker'}});


profile viewer