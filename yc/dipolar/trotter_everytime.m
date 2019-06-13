%20190508
% compute function ZZ(t) and XX(t) of dipolar_z with pi decoupling,
% consider pulse width
% at T = \infty,
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/8

%% runtime parameters
N_atom = 10;

alpha = 3;
hx = 1;
Jz = 1;
u = 1;
hz = 0;

Jza = u*Jz; % * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [1, -1/2, -1/2];
bc = 'p';
sym = ['kp'];

dtau = 0.5;
tau_list = [dtau:dtau:1.5];

t_end = 5e1;

load_save = false;
filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau05.mat',N_atom, hx, hz);
     
profile on
profile clear

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'x',1/2);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'z',1/2);
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
opers = {OperatorClass(N_atom,'x',1/2), OperatorClass(N_atom,'z',1/2), OperatorClass(N_atom,'y',1/2)};
for oper_index = 1: length(opers)
    opers{oper_index}.symmetrize(sym);
end

%% begin of for
[ZZ_list,ZZ_ideal, ZZ_dev] = deal( cell(1, length(tau_list)+1) );

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    Ut = H2U(Hz, tau)* H2U(Hx, tau);
    %Ut = H2U(Hz* tau, 1) * H2U(pi* Hx + Hz* 0., 1);
    period_N = floor(t_end/tau);
    
    [Ut_t, Ut_tau0] = deal( H2U(Hz, 0));
    Utau0 = H2U(Hz+Hx, tau);
    [ZZ_list{tau_index}, ZZ_dev{tau_index}] = deal( zeros(length(opers), period_N) );
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

figure(43);
clf;

subplot(1,2,1);
hold on
legends{end+1} = ['$J_z\tau=0$'];
ycPlot([1: t_end], ZZ_ideal{2}(3,:), '', '', '' ,[]);

for t_index = 1:length(tau_list)
    legends{end+1} = ['$J_z\tau=',mat2str( tau_list(t_index) ),'$'];
end
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_list{tau_index},2)], ZZ_list{tau_index}(3,:) ...
        , '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['$Dip_z+Y, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
ylim([-0.5 0.5])

subplot(1,2,2);
hold on
legends{end+1} = ['$J_z\tau=0$'];
ycPlot([1: t_end], ZZ_ideal{2}(1,:), '', '', '' ,[]);

for t_index = 1:length(tau_list)
    legends{end+1} = ['$J_z\tau=',mat2str( tau_list(t_index) ),'$'];
end
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_list{tau_index},2)], ZZ_list{tau_index}(1,:) ...
        , '$J_zt$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
        ['$Dip_z+Y,\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
ylim([-0.2 0.2])

%{
for t_index = 1:length(tau_list)
    legends{end+1} = ['$\Delta ZZ, J_z\tau=',mat2str( tau_list(t_index) ),'$'];
end
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_dev{tau_index},2)], ZZ_dev{tau_index}(1,:) ...
        , '$J_zt$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
%
%ylim([-0.2,0.2])


subplot(1,2,2)
hold on
%ycPlot([1: t_end], ZZ_list{end}(2,:), '', '', '' );
for tau_index = 1: length(tau_list)
    ycPlot(tau_list(tau_index)* [1:size(ZZ_list{tau_index},2)], ZZ_list{tau_index}(3,:) ...
        , '$J_zt$', '$<X(t)X(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
end
%}

profile viewer









