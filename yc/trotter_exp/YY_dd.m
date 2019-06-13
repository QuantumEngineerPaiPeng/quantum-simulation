%20190429
% try to produce the dipolar longtime XX plot from experiment data with DD
% Created by Chao Yin
% last update: 4/29

%% runtime parameters
exp_num_list = [40:49, 51];
exp_name = 'expt3';

load_data = true;

mean_N = 10;

%% load data
if load_data
    YY = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY{end+1} = -real(B);
    end
end
%% fixt
tau_N = length(exp_num_list);
htau_list = [15:15:165];

fixt_n_list = [];
fixt_YY = [];
for tau_index = 2:tau_N
    fixt_n_list(end+1) = length(YY{tau_index});
    fixt_YY(end+1) = mean( YY{tau_index}(end-4:end) );
end

%% plot
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
        ['$h=J, \pi_X$-pulse DD'], legends);
end

subplot(1,2,2)
hold on
legends =  {'experiment','numeric, $h\tau$'};
ycPlot(htau_list(2:end), fixt_YY, '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends, {{'Marker'}});
ycPlot(tau_list *180/pi , 2*fixt_YY(1)/(XX_list(end,2)+ XX_list(end,3)) * XX_list(end,:), '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 30^\circ$'], legends);

%% plot all YY, varying tau and t
figure(1)
clf;

subplot(1,2,1);
legends = cell(0);
legends{end+1} = 'fix $t$';
for tau_index = fixt_start: fixt_end
    legends{end+1} = ['$h\tau=',mat2str( htau_list(tau_index) ),'^{\circ}$'];
end

scatter(fixt_n_list, fixt_YY, 'LineWidth', 2);
hold on
ycPlot(n_list, YY(:,fixt_start: fixt_end), '$n_{period}$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, pulse=\pi+h\tau$'], legends);


subplot(1,2,2)
hold on
legends =  {'experiment','numeric'};
ycPlot(fixt_htau, fixt_YY, '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$'], legends, {{'Marker'}});
ycPlot(-tau_list *180/pi , fixt_YY(end)/XX_list(end,1) * XX_list(end,:), '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 15^\circ$'], legends, {{'Marker'}});
%xlim([0 180])

%% plot YY vs tau, fix n=n_max & YY vs tau, fix t=tau_min*n_max
figure(2)
clf;
subplot(1,2,1)
ycPlot(htau_list, mean(YY(end-mean_N: end, :), 1), '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$']);

subplot(1,2,2)
hold on
ycPlot(fixt_htau, fixt_YY, '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$']);
ycPlot(tau_list, fixt_YY(1)* XX_list, '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$'], {'experiment','numeric'});

profile viewer