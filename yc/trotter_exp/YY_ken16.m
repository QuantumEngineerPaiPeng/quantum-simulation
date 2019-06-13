%20190514
% try to produce the dipolar longtime XX plot from experiment data with Ken
% pulse
% Created by Chao Yin
% last update: 5/14

%% runtime parameters
exp_num_list = [149:159];
exp_name = 'expt3';

load_data = true;

mean_N = 10;

%% load data
if load_data
    YY = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY{end+1} = real(B);
    end
end

%% fixt
tau_N = length(exp_num_list);
htau_list = [15:15:165];

fixt_n_list = [];
fixt_YY = [];
for tau_index = 2:tau_N
    fixt_n_list(end+1) = length(YY{tau_index});
    fixt_YY(end+1) = mean( YY{tau_index}(end:end) );
end

%% plot
exp_name = 'Ken 16-pulse';
figure(5)
clf;
subplot(1,2,1)
legends = {'fix $t$'};
for tau_index = 1:tau_N
    legends{end+1} = ['$h\tau=',mat2str( htau_list(tau_index) ),'^{\circ}$'];
end

hold on
scatter(fixt_n_list, fixt_YY, 'LineWidth', 2);
for tau_index = 1: tau_N 
    ycPlot([1: length(YY{tau_index})], YY{tau_index}, '$n_{period}$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
        ['$h=J$, ',exp_name], legends);
end

subplot(1,2,2)
hold on
legends =  {'experiment','numeric, $h\tau$'};
ycPlot(htau_list(2:end), fixt_YY./YY{1}(1), '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(tau_list *180/pi , ...%2*fixt_YY(1)/(XX_list(2)+ XX_list(3)) * 
    XX_list, '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});

%% all good experiments
legends = {'no dd, $\tau_1=2.95$', 'no dd, $\tau_1=1.5$'...
    ,'$\pi$ dd', 'Ken 16-pulse', '$\pi$ dd with echo','numeric $N=16$'}

figure(90)
clf;

normal_factor = (XX_list(2)+XX_list(3) )/2;
ycPlot(htau_list_nopi, fixt_YY_nopi{1}./fixt_YY_nopi{1}(1) *XX_list(1), '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(htau_list_nopi, fixt_YY_nopi{2}./fixt_YY_nopi{2}(1) *XX_list(1), '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(htau_list_dd, fixt_YY_dd./fixt_YY_dd(1) *normal_factor, '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(htau_list_ken, fixt_YY_ken./fixt_YY_ken(1) *normal_factor, '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
normal_factor = XX_list(1);
ycPlot(htau_list_dd_echo, fixt_YY_dd_echo./fixt_YY_dd_echo(1) *normal_factor, '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(tau_list *180/pi , ...%2*fixt_YY(1)/(XX_list(2)+ XX_list(3)) * 
    XX_list(end, :), '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J$'], legends, { {'color', 'k'}, {'line_type', ':'}});

ylim([0 0.9])




