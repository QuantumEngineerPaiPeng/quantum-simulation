%20190331
% compute partition function Z of 2d ising model
% Z = \sum e^{-H}
% H = K1 \sigma_j,k \sigma_j+1,k + K2 \sigma_j,k \sigma_j,k+1 + h
% \sigma_j,k
% Created by Chao Yin
% last change: 3/31

%% runtime parameters
N1 = 3;
N2 = 3;

hx = 1;
Jz = 1;
hz = 0;

dtau = 0.02;
tau_list = [dtau: dtau: 1];
%method = {'monte', 10^6};
method = {'enumerate'};

load_save = false;
exact = true;

profile on
profile clear

%% real parameters, solvable
cri_htau = 1/2* log(1+sqrt(2));
Z_list = tau_list;
for tau_index = 1: length(tau_list)
    tau = tau_list(tau_index);
    Z_list(tau_index) = par_fun_ising( N1, N2, -tau, -tau, 0, method );
end % for tau_index = 1: length(tau_list)
free_e = log(2) + log( Z_list) ./ (N1* N2);

%% exact free energy
if exact
    cri_htau = 1/2* log(1+sqrt(2));
    free_exact = tau_list;
    for tau_index = 1: length(tau_list)
        tau = tau_list(tau_index);
        func = @(x,y) log( cosh(2*tau) * cosh(2* tau)- ...
            sinh(2*tau)* cos(x) - sinh(2*tau)* cos(y) );
        free_exact(tau_index) = log(2) + 1/8/pi^2 * integral2(func, -pi, pi, -pi, pi);
    end % for tau_index = 1: length(tau_list)
end % if exact

%% solvable ising plot
legend = {['$N=',mat2str(N1), '\times',mat2str(N2),'$',', exact'], '$N=\infty$, analytical', 'transition point'};
extra = { {'x_par', cri_htau} };

figure(52);
clf;
%' ,N_{sample}=',sprintf('%0.5g',method{2}),...
subplot(1,2,1)
ycPlot(tau_list, [free_e; real( free_exact)], '$K_1$', '$\ln{Z}/N$', ...,
    ['ising, $h=', mat2str(0),', K_1+2K_2=3K_c' ,'$'], legend, extra);

subplot(1,2,2)
ycPlot(tau_list(2:end-1), [diff(free_e, 2); diff(free_exact, 2)], '$K_1$', '$ \partial^2_{K_1} \ln{Z}/N $', ...,
    ['ising, $h=', mat2str(0),', K_1+2K_2=3K_c' ,'$'], legend, extra);

profile viewer



