%20190426
% try to produce the dipolar longtime XX plot from experiment data
% Created by Chao Yin

%% runtime parameters
exp_num_list = [204 205]; % XX, ZZ
exp_name = 'expt3';

load_data = true;
run_simu = true;

%% load data
if load_data
    YY = cell(0);
    for exp_index = 1: length(exp_num_list)
        data_analysis(exp_num_list(exp_index), exp_name);
        YY{end+1} = real(B);
    end
end
n_N = size(YY{1}, 1);
n_list = [1: n_N];
tau_N = size(YY{2}, 2);
htau_list = [15:15:165]; % [ [-165:15:-15], [15:15:90] ];


%% fluctuation
XX_mean = mean(YY{1}(end/2:end, :), 1);
ZZ_mean = mean(YY{2}(end/2:end, :), 1);
XX_var = var(YY{1}(end/2:end, :), 1);
ZZ_var = var(YY{2}(end/2:end, :), 1);

%% generate ideal numerics
if run_simu
N_atom = 8;
alpha = 3;
hx = 1;
Jz = 1;
Jza = Jz;
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau = 15;
tau_list = [dtau:dtau:165];
dtau = dtau*pi/180;
tau_list = tau_list*pi/180;

n_end = n_N*2;

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'y',1/2);
Hx.symmetrize(sym);

%% test operator 
opers = {OperatorClass(N_atom,'x',1/2), OperatorClass(N_atom,'y',1/2), OperatorClass(N_atom,'z',1/2)};
for oper_index = 1: length(opers)
    opers{oper_index}.symmetrize(sym);
end

%% begin of for
[ZZ_list, ZZ_ideal, ZZ_dev] = deal( cell(1, length(tau_list)) );

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    Ut = H2U(Hz, tau)* H2U(Hx, tau);
    %Ut = H2U(Hz* tau, 1) * H2U(pi* Hx + Hz* 0., 1);
    period_N = n_end;
    
    [Ut_t, Ut_tau0] = deal( H2U(Hz, 0));
    Utau0 = H2U(Hz+Hx, tau);
    [ZZ_list{tau_index}, ZZ_dev{tau_index}, ZZ_ideal{tau_index}] = deal( zeros(length(opers), period_N) );
    for t_index = 1: period_N
        Ut_t = Ut* Ut_t;
        Ut_tau0 = Utau0* Ut_tau0;
        for oper_index = 1: length(opers)
            ZZ_list{tau_index}(oper_index, t_index) = trace(Ut_t' * opers{oper_index} * Ut_t * opers{oper_index});
            ZZ_ideal{tau_index}(oper_index, t_index) = trace(Ut_tau0' * opers{oper_index} * Ut_tau0 * opers{oper_index});
        end
    end % for t_index = 1: length(t_list)
    ZZ_list{tau_index} = real( ZZ_list{tau_index})/ 2^N_atom *4/ N_atom;
    ZZ_ideal{tau_index} = real( ZZ_ideal{tau_index})/ 2^N_atom *4/ N_atom;
    ZZ_dev{tau_index} = ZZ_list{tau_index} - ZZ_ideal{tau_index};

end % for tau_index = 1:length(tau_list)
end % if run_simu

%% raw data
figure(13)
clf;

hold on
legends = {};
for tau_index =1:tau_N
    legends{end+1} = ['$\tau=',mat2str(htau_list(tau_index)),'$'];
    compare_dots = [tau_index:tau_index:n_N];
    lines = ycPlot(compare_dots, YY{1}(1:length(compare_dots), tau_index), ...
        '$Jt/15^\circ$', '$<X(n\tau)X(0)>_{\beta=0}$', ...,
    ['$h=J, pulse=h\tau$'],legends);
    colors = get(lines, 'color');
    %ycPlot(n_list, YY{2}(:, tau_index) - YY{2}(:,1), '$n_{period}/2$', '$\Delta<S(n\tau)S(0)>_{\beta=0}$', ...,
       % ['$h=J, pulse=h\tau, ideal: \tau=15$'], legends, {{'line_type','--'}, {'color', colors}});
end
%ylim([-2000 1000])

%% all YY and fix tfor all exps in one subplot 
figure(11)
clf;

subplot(1,2,1)
hold on
legends = {};
deviation = [];
for tau_index = 1:tau_N
    legends{end+1} = ['$\tau=',mat2str(htau_list(tau_index)),'$'];
    compare_dots = [tau_index:tau_index:n_N];
    temp = (YY{1}(1:length(compare_dots), tau_index) - YY{1}(compare_dots,1))./ abs( YY{1}(compare_dots,1) );
    deviation(end+1) = sqrt( mean( temp.^2));
    lines = ycPlot(compare_dots, temp, ...
        '$n_{period}/2$', '$\Delta<X(n\tau)X(0)>_{\beta=0}$', ...,
    ['$h=J, pulse=h\tau$'],legends);
    colors = get(lines, 'color');
    %ycPlot(n_list, YY{2}(:, tau_index) - YY{2}(:,1), '$n_{period}/2$', '$\Delta<S(n\tau)S(0)>_{\beta=0}$', ...,
       % ['$h=J, pulse=h\tau, ideal: \tau=15$'], legends, {{'line_type','--'}, {'color', colors}});
end
%xlim([0 30])
%ylim([-2000 1000])

% assume tau=15 is ideal

subplot(1,2,2)
hold on
legends =  {'XX'};
ycPlot(htau_list,deviation ...%fixt_YY{1}(1) *XX_list(1)
    , '$h\tau/^\circ$', '', ...,
    ['relative $\sqrt{ \langle \Delta<X(n\tau)X(0)>_{\beta=0}^2} \rangle_n$'], legends, {{'Marker'}});
%set(gca, 'yscale', 'log')
%}

%% XX deviation from ideal numerics
figure(12)
clf;

subplot(1,2,1)
hold on
legends = {'$XX^N, \tau=0$', '$XX^N, \tau=15$', '$XX^E, \tau=15$',...
    '$XX^N, \tau=0$', '$XX^N, \tau=90$', '$XX^E, \tau=90$'};
for tau_index = 1:6:7
    lines = ycPlot(n_list, ZZ_ideal{tau_index}(1,2:2:end) ...
        , '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
    colors = get(lines, 'color');
    ycPlot(n_list, ZZ_list{tau_index}(1,2:2:end), '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends, {{'line_type',':'}, {'color', colors}});
    ycPlot(n_list, YY{1}(:, tau_index)./ YY{1}(1, tau_index)* ZZ_list{tau_index}(1,2),...
        '$n_{period}/2$', '$<X(n\tau)X(0)>_{\beta=0}$', ...,
    ['XX, $h=J, pulse=h\tau$'],legends ,{{'line_type','--'}, {'color', colors}});
end
%ylim([-1000 3000])

subplot(1,2,2)
hold on
legends =  {'XX','ZZ'};
ycPlot(htau_list, [ sqrt(XX_var); sqrt(ZZ_var)] ...%fixt_YY{1}(1) *XX_list(1)
    , '$h\tau$', '$\Delta SS$', ...,
    ['$h=J, n_{period}=64:128$'], legends, {{'Marker'}});

profile viewer