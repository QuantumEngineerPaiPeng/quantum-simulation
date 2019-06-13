%20190507
% try to produce the dipolar longtime XX plot from experiment data
% Created by Chao Yin

%% runtime parameters
exp_num_list = [34];
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
YY = YY{1};
n_N = size(YY, 1);
n_list = [1: n_N];
tau_N = size(YY, 2);
htau_list = [ [-165:15:-15], [15:15:90] ];


%% fix t data
[fixt_YY, fixt_n_list] = deal( [] );
fixt_htau = htau_list;
fixt_n_list = round( [ n_N ./ [11:-1:1], n_N ./ [1:6] ] );
fixt_YY = zeros(1, length(fixt_htau));
for fixt_index = 1: length(fixt_htau)
    fixt_YY(fixt_index) = YY( round(fixt_n_list(fixt_index)), fixt_index);
end

%% all YY and fix-t for all exps in one subplot 
figure(12)
clf;
subplot(1,2,1)
legends = {'fix $t$'};
for tau_index = 1:tau_N
    %legends{end+1} = ['$h\tau=',mat2str( htau_list(tau_index) ),'^{\circ}$'];
end

hold on
scatter(fixt_n_list, fixt_YY, 'LineWidth', 2);
for tau_index = 1: tau_N 
    ycPlot(n_list, YY, '$n_{period}/2$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
        ['$h=J, pulse=\pi+h\tau$'], legends);
end

subplot(1,2,2)
hold on
legends =  {'experiment','numeric'};
ycPlot(htau_list, [fixt_YY(1:11)./fixt_YY(11) *XX_list(1), fixt_YY(12:end)./fixt_YY(12) *XX_list(1) ], '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 15^\circ$'], legends, {{'Marker'}});
plt_x = [ -tau_list(end:-1:1) *180/pi, tau_list *180/pi] ;
ycPlot( plt_x , [XX_list(end:-1:1), XX_list], '$h\tau$', '$<Y(n\tau)Y(0)>_{\beta=0}$', ...,
    ['$h=J, ht=128\times 15^\circ$'], legends);

profile viewer


