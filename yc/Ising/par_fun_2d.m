%20190331
% compute partition function Z of 2d ising model
% Z = \sum e^{-H}
% H = K1 \sigma_j,k \sigma_j+1,k + K2 \sigma_j,k \sigma_j,k+1 + h
% \sigma_j,k
% in our case, K_1 = -i\tau Jz, K_2= i\pi/4 - 1/2 ln(tan(hx \tau)), 
% h = -i\tau hz
% Created by Chao Yin
% last change: 3/31

%% runtime parameters
N1 = 4;
N2 = 4;

hx = 1;
Jz = 1;
hz = 0;

dtau = 0.02;
tau_list = [2: dtau: 3];
%method = {'monte', 10^6};
method = {'enumerate'};

load_save = false;
exact = true;

profile on
profile clear


%% calculate Z
Z_list = tau_list;
for tau_index = 1: length(tau_list)
    tau = tau_list(tau_index);
    K1 = -i* tau* Jz;
    K2 = i*pi/4 - 1/2* log(tan(hx* tau));
    h = -i* tau* hz;
    Z_list(tau_index) = par_fun_ising( N1, N2, K1, K2, h, method );
end % for tau_index = 1: length(tau_list)
free_e = log(2) + log( Z_list) ./ (N1* N2);

%% exact free energy
if exact
    free_exact = tau_list;
    for tau_index = 1: length(tau_list)
        tau = tau_list(tau_index);
        K1 = -i* tau* Jz;
        K2 = i*pi/4 - 1/2* log(tan(hx* tau));
        func = @(x,y) log( cosh(2*K1)* cosh(2*K2) - sinh(-2*K1)*cos(x) - sinh(-2*K2)* cos(y) );
        free_exact(tau_index) = log(2) + 1/8/pi^2 * integral2(func, -pi, pi, -pi, pi);
    end % for tau_index = 1: length(tau_list)
end % if exact

%% exact result of ising at hx=0
cri_htau = pi ./ (1+ 1./ (2* hx) );

%% solvable ising plot
legend = {['$N=',mat2str(N1), '\times',mat2str(N2),'$',', exact'], '$N=\infty$, analytical', 'transition point'};
extra = { {'x_par', cri_htau} };

figure(51);
clf;
subplot(1,2,1)
ycPlot(tau_list, [free_e; free_exact], '$h_x\tau$', '$\ln{Z}/N$', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx), ' ,\frac{h_z}{J_z}=',mat2str(hz) , '$'], legend, extra);

subplot(1,2,2)
ycPlot(tau_list(2:end-1), [diff( real(free_exact), 2)], '$h_x\tau$', '$ C_v $', ...,
    ['ising,$\frac{h_x}{J_z}=',mat2str(hx), ' ,\frac{h_z}{J_z}=',mat2str(hz) , '$'], legend, extra);

profile viewer



