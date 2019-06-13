%20190509
% random guess to solve the problem in PHYS. REV. X 8, 031086 (2018)
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/9

%% runtime parameters
N_atom = 10;

g = 1;
hx = [-2 2];

ha = 4;

model = [0, 0, 1];
bc = 'p';
cplist = [1];
sym = [];
find_state = false;

T = 3;
dt = 0.05;
Nt = floor(T/ dt);

nEpi = 1e3;
n_run = 5;

profile on
profile clear

%% find the initial and final state as the GS of traverse ising
H0 = OperatorClass(N_atom,model,-1/4,bc,cplist) - g* OperatorClass(N_atom,'z',1/2) ;
Hx = OperatorClass(N_atom,'x',1/2);

H_ini = H0 - hx(1) * Hx;
H_fin = H0 - hx(2) * Hx;

if find_state
H_ini.diagonalize();
psi_ini = H_ini.eigsys{1}.V(:, 1);
max(abs( H_ini.matrix{1}* psi_ini - H_ini.eigsys{1}.D(1)* psi_ini ))

H_fin.diagonalize();
psi_fin = H_fin.eigsys{1}.V(:, 1);
max(abs( H_fin.matrix{1}* psi_fin - H_fin.eigsys{1}.D(1)* psi_fin ))
end

%%
Ut = { H2U( H0 - ha*Hx, dt), H2U( H0 + ha*Hx, dt) };
F_list = zeros(n_run, nEpi);

for run_index = 1: n_run
    for epi_index = 1: nEpi
        psi = psi_ini;
        pulse = randi(2, 1, Nt);
        for t_index = 1: Nt
            psi = Ut{pulse(t_index)}.matrix{1} * psi;
        end
        F_list(run_index, epi_index) = abs( psi'* psi_fin )^2;
    end % for epi_index = 1: nEpi
end % for run_index = 1: n_run

%%
set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultColorbarTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextFontSize',17)
set(groot,'defaultAxesFontSize',17)

figure(51)
clf;
hold on

for run_index = 1: n_run
    scatter( [1:nEpi], F_list(run_index, :), 'LineWidth', 2);
end
x_lim = get(gca,'Xlim');
plot( x_lim, 0.9*[1 1],'-','Linewidth',2 );
xlabel('episode')
ylabel('Fidelity')
title('PRX random guess, 5 colors: 5 runs')

profile viewer


