%20190426
% try to produce the dipolar longtime XX plot from experiment data
% Created by Chao Yin

%% runtime parameters
exp_num_list = [33,50]; % 191,194
exp_name = 'expt3';

load_data = true;

mean_N = 10;

%% load data
if load_data
    YY = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY{end+1} = imag(B);
    end
end
n_N = size(YY{1}, 1);
n_list = [1: n_N];
tau_N = size(YY{2}, 2);
htau_list = {[15:15:150],[15:15:150]}; % [ [-165:15:-15], [15:15:90] ];


%% fix t data
[fixt_YY, fixt_n_list] = deal( cell(1, length(exp_num_list)) );
for exp_index = 1: length(exp_num_list)
    fixt_htau = htau_list{exp_index};
    fixt_n_list{exp_index} = round( n_N ./ [1: length(fixt_htau)] );
    fixt_YY{exp_index} = zeros(1, length(fixt_htau));
    for fixt_index = 1: length(fixt_htau)
        fixt_YY{exp_index}(fixt_index) = YY{exp_index}( round(fixt_n_list{exp_index}(fixt_index)), fixt_index);
    end
end % for exp_index = 1: length(exp_num_list)
%% all YY and fix tfor all exps in one subplot 
figure(11)
clf;

subplot(1,2,1)
hold on
legends = {};
legends{end+1} = 'fix $t$, $\tau_1=2.95\mu s$';
legends{end+1} = 'fix $t$, $\tau_1=1.5\mu s$';
for tau_index = 1:tau_N
    %legends{end+1} = ['$h\tau=',mat2str( htau_list{2}(tau_index) ),'^{\circ}$'];
end
scatter(fixt_n_list{1}, fixt_YY{1}, 'LineWidth', 2);
scatter(fixt_n_list{2}, fixt_YY{2},'^', 'LineWidth', 2);
lines = ycPlot(n_list, YY{1}, '$n_{period}/2$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, pulse=h\tau$'],legends);
colors = get(lines, 'color');
for tau_index = 1:tau_N
    ycPlot(n_list, YY{2}(:, tau_index), '$n_{period}/2$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
        ['$h=J, pulse=h\tau$'], legends, {{'line_type','--'}, {'color', colors{tau_index}}});
end


subplot(1,2,2)
hold on
legends =  {'expt, $\tau_1=2.95\mu s$', 'expt, $\tau_1=1.5\mu s$','numeric'};
ycPlot(fixt_htau, fixt_YY{1}./ YY{1}(1,1) ...%fixt_YY{1}(1) *XX_list(1)
    , '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$'], legends, {{'Marker'}});
ycPlot(fixt_htau, fixt_YY{2}./ YY{2}(1,1) ...%fixt_YY{1}(1) *XX_list(1)
    , '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, n_{period}=128$'], legends, {{'Marker'}});
ycPlot(tau_list *180/pi , XX_list, '$h\tau/^\circ$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 15^\circ$'], legends);

profile viewer