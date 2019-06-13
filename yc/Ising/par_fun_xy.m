%20190401
% compute partition function Z of 2d ising model
% Z = \sum e^{-H}
% H = K1 \sigma_j,k \sigma_j+1,k + K2 \sigma_j,k \sigma_j,k+1 + h
% \sigma_j,k
% in our case, K_1 = -i x, K_2= i\pi/4 + y, 
% h = 0
% Created by Chao Yin
% last change: 4/01

%% runtime parameters
N1 = 4;
N2 = 4;

hx = 1;
Jz = 1;
hz = 0;

x2 = pi/3; % x2 = 2*x
dy = 0.01;
y_list = [-0.5: dy: 0];
%method = {'monte', 10^6};
method = {'enumerate'};

load_save = false;
exact = true;

profile on
profile clear


%% calculate Z
Z_list = y_list;
for y_index = 1: length(y_list)
    y = y_list(y_index);
    K1 = -i* x2/2;
    K2 = i*pi/4 +y;
    h = 0;
    Z_list(y_index) = par_fun_ising( N1, N2, K1, K2, h, method );
end % for tau_index = 1: length(tau_list)
free_e = log(2) + log( Z_list) ./ (N1* N2);

%% exact free energy
if exact
    free_exact = y_list;
    for y_index = 1: length(y_list)
        y = y_list(y_index);
        K1 = -i* x2/2;
        K2 = i*pi/4 +y;
        func = @(theta1,theta2) log( cosh(2*K1) * cosh(2* K2)- ...
            sinh(2*K1)* cos(theta1) - sinh(2*K2)* cos(theta2) );
        %func = @(theta1, theta2) log( cos(x2)* sinh(2*y) - cosh(2*y)*cos(theta1) + sin(x2)* cos(theta2) );
        free_exact(y_index) = log(2) + 1/8/pi^2 * integral2(func, -pi, pi, -pi, pi);
    end % for tau_index = 1: length(tau_list)
end % if exact

%% exact result of ising at hx=0
cri_y = -1/2* acosh(1/ sin(x2));

%% plot
legend = {['$N=',mat2str(N1), '\times',mat2str(N2),'$',', exact'], '$N=\infty$, analytical', 'transition point'};
extra = { {'x_par', cri_y} };

figure(53);
clf;
subplot(1,2,1)
ycPlot(y_list, [free_e; free_exact], '$y$', '$\ln{Z}/N$', ...,
    ['Ising, $K_1=-i\pi/3, K_2=i\pi/4 + y$'], legend, extra);

subplot(1,2,2)
ycPlot(y_list(2:end), [diff( real(free_e), 1); diff( real(free_exact), 1)], '$y$', '$ \partial_{y} \ln{Z}/N $', ...,
    ['Ising, $K_1=-i\pi/3, K_2=i\pi/4 + y$'], legend, extra);

profile viewer



