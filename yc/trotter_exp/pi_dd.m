%20190501
% try to produce the dipolar longtime XX plot from experiment data with DD
% Created by Chao Yin
% last update: 5/01

%% runtime parameters

exp_num_list = [  [72:83]];
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
    
exp_name = 'expt3';

load_data = true;

%% load data
if load_data
    YY = cell(0);
    YX = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY{end+1} = real(B);
        YX{end+1} = imag(B);
    end
end
%% fixt
tau_N = length(exp_num_list);
n_N = length(YY{1});
J = 33e-3;

tau_list = (d4+ p2);
Jtau_list = J* tau_list;

fixt_n_list = round( tau_list(1)*n_N ./ tau_list );
fixt_YY = [];
for tau_index = 1:tau_N
    fixt_YY(end+1) = YY{tau_index}(fixt_n_list(tau_index));
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