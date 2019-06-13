%20190318
% compute observable ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/18

%% runtime parameters
N_atom = 4;

alpha = 3;
hx = 1;
Jz = 1;
u = 1;
hz = 0;

Jza = u*Jz; % * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau = 0.2;
tau_list = [dtau:dtau:4];

t_list = [10:10:40];

load_save = false;
filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau02.mat',N_atom, hx, hz);
     
profile on
profile clear

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
Z0 = OperatorClass(N_atom,'z',1/2);
X0 = OperatorClass(N_atom,'x',1/2);
Z0.symmetrize(sym);
X0.symmetrize(sym);

%% begin of for
ZZ_list = zeros(length(t_list)+1, length(tau_list));
XX_list = zeros(length(t_list)+1, length(tau_list));

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    if loaded
        Ut = Ut_all{tau_index};
    else
        Ut = U1 * Ut * U2;
        Ut_all{end+1} = copy( Ut);
    end
    
    for t_index = 1: length(t_list)+1
        if t_index == length(t_list)+1
            t_interval = 'ED';
        else
            t_interval = [t_list(t_index), t_list(t_index)];
        end
    
        temp =  ZZ_longtime(Ut_all{tau_index}, t_interval, {Z0,X0});
        ZZ_list(t_index, tau_index) = temp(1, :)*4/N_atom;
        XX_list(t_index, tau_index) = temp(2, :)*4/N_atom;
    
    end % for t_index = 1: length(t_list)

end % for tau_index = 1:length(tau_list)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end

%% plot ZZ, XX vs N_atom
legends = cell(0);
for t_index = 1:length(t_list)
    legends{end+1} = ['$J_zt=',mat2str( t_list(t_index) ),'$'];
end
legends{end+1} = ['$J_zt=\infty$'];

figure(41);
clf;

subplot(1,2,1);
ycPlot(tau_list, ZZ_list, '$J_z\tau$', '$<Z(t)Z(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);
hold on
subplot(1,2,2);
ycPlot(tau_list, XX_list, '$J_z\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);

%% only XX (hz=0)
figure(42)
ycPlot(tau_list, XX_list, '$J_z\tau$', '$<X(t)X(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),'$'], legends);


profile viewer









