%20190508
% compute function ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty, and its fluctuation, deviation
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5.8

%% runtime parameters
N_atom = 10;

alpha = 3;
hx = 1;
Jz = 1;
u = 1;
hz = 0;
pi_w = 0.;
angle_e = 10; % degree
phase_e = 0;

Jza = u*Jz; % * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau =7.6/180*pi;
tau_list = [dtau:dtau:1.2];

t_end = 5e1;

load_save = false;
filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau05.mat',N_atom, hx, hz);
     
profile on
profile clear

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);

Hx=hx* (cos(phase_e*pi/180)* OperatorClass(N_atom,'x',1/2) + sin(phase_e*pi/180) ...
    * OperatorClass(N_atom,'y',1/2));
Hx.symmetrize(sym);

if load_save && exist(filename,'file')
    load(filename);
    loaded = true
else
    Ut_all = cell(0);
    Ut=OperatorClass(N_atom);
    Ut.matrix = { speye(2^N_atom) };
    Ut.symmetrize(sym);
    U1 = H2U(Hz, dtau);
    U2 = H2U(Hx, dtau);
    
    loaded = false
end % if exist(filename,'file')

%% test operator 
opers = {OperatorClass(N_atom,'z',1/2), OperatorClass(N_atom,'x',1/2)};
for oper_index = 1: length(opers)
    opers{oper_index}.symmetrize(sym);
end

%% begin of for
[ZZ_list, ZZ_dev] = deal( cell(1, length(tau_list)+1) );

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    Ut = H2U(Hz,tau) * H2U((pi+ angle_e*pi/180)* Hx + Hz* pi_w,1);
    period_N = floor(t_end/tau);
    
    [Ut_t, Ut_tau0] = deal( H2U(Hz, 0));
    Utau0 = H2U(Hz+Hx, tau);
    [ZZ_list{tau_index}, ZZ_dev{tau_index}] = deal( zeros(length(opers), period_N) );
    for t_index = 1: period_N
        Ut_t = Ut* Ut_t;
        Ut_tau0 = Utau0* Ut_tau0;
        for oper_index = 1: length(opers)
            ZZ_list{tau_index}(oper_index, t_index) = trace(Ut_t' * opers{oper_index} * Ut_t * opers{oper_index});
            ZZ_dev{tau_index}(oper_index, t_index) = ZZ_list{tau_index}(oper_index, t_index) ...
                - trace(Ut_tau0' * opers{oper_index} * Ut_tau0 * opers{oper_index});
        end
    end % for t_index = 1: length(t_list)
    ZZ_list{tau_index} = real( ZZ_list{tau_index})/ 2^N_atom *4/ N_atom;
    ZZ_dev{tau_index} = real( ZZ_dev{tau_index})/ 2^N_atom *4/ N_atom;

end % for tau_index = 1:length(tau_list)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end

%% deviation from tau=0 and its variance, then plot
ZZ_dev_var = zeros(length(opers), length(tau_list));
for tau_index = 1: length(tau_list)
    ZZ_dev_var(:, tau_index) = var(ZZ_dev{tau_index}, 0, 2);
end
figure(45)
ycPlot(tau_list, ZZ_dev_var(1,:) ...
        , '$J_z\tau$', '$\Delta<Z(t)Z(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], []);

%% plot ZZ, XX vs t
legends = cell(0);

figure(44);
clf;
%subplot(1,2,1);
hold on

%ycPlot([1: t_end], ZZ_list{end-1}(1,:), '', '', '' ,[]);
%tau_list = [15, 30, 45, 60]
tau_list_angle = [15, 30, 45, 60];
for t_index = 1:length(tau_list)
%    legends{end+1} = ['$J_z\tau=',mat2str( tau_list_angle(t_index) ),'^\circ$'];
end
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_list{tau_index},2)], ZZ_list{tau_index}(2,:) ...
        , '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['$J\tau_\pi=', mat2str(pi_w) ...
        , ' ,\Delta angle=',mat2str(angle_e), '^\circ ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
ylim([-0.2 1])
%}

%{
subplot(1,2,2)
hold on
ycPlot([1: t_end], ZZ_list{end}(2,:), '', '', '' );
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_list{tau_index},2)], ZZ_list{tau_index}(2,:) ...
        , '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
%}

profile viewer









