%20190319
% compute time fluctuation & deviation of ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty,
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/19

%% runtime parameters
N_all = [6:2:12];

alpha = 3;
hx = 1;
Jz = 1;
u = 1;
hz = 0;

Jza = u*Jz; % * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau = 0.1;
tau_list = [dtau:dtau:3];

t_end = 5e1;

load_save = false;
     
profile on
profile clear

%% begin of diff N
ZZ_dev_var = zeros(length(N_all), length(tau_list));
for N_index= 1: length(N_all)
N_atom = N_all(N_index);
cplist = Jza * [1:N_atom-1].^(-alpha) ;
filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau02.mat',N_atom, hx, hz);

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'x',1/2);
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
opers = {OperatorClass(N_atom,'z',1/2)};
for oper_index = 1: length(opers)
    opers{oper_index}.symmetrize(sym);
end

%% begin of for
[ZZ_list, ZZ_dev] = deal( cell(1, length(tau_list)+1) );

for tau_index = 1:length(tau_list)+1
    if tau_index == length(tau_list)+1
        tau = 1;
        Ut = H2U(Hz+Hx, 1);
    else
        tau = tau_list(tau_index);
        if loaded
            Ut = Ut_all{tau_index};
        else
            Ut = U1 * Ut * U2;
            Ut_all{end+1} = copy( Ut);
        end

    end % if tau_index == 0
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
for tau_index = 1: length(tau_list)
    ZZ_dev_var(N_index, tau_index) = sqrt( mean(ZZ_dev{tau_index}.^2) );
    %mean(ZZ_dev{tau_index})
end

end % for N_index


%% plot 
legends = cell(0);
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end

figure(46)
ycPlot(tau_list, ZZ_dev_var ...
        , '$J_z\tau$', '$\Delta<Z(t)Z(0)>_{\beta=0}$', ...,
        ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,t_{end}=', mat2str(t_end),'$'], legends);
set(gca,'yscale','log')

profile viewer









