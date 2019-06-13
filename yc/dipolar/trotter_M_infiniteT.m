%20190315
% compute observable ZZ(t->\infty) and XX of dipolar_z with hx and hz
% interacions at T = \infty
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N, Z_0 = S^z / N
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 8;

alpha = 3;
hx = 1;
Jz = 1;
u = 1;
hz = 1;

Jza = u*Jz; % * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau = 0.1;
tau_list = [dtau:dtau:6];

F_N = 10^2;

load_save = false;
filename=sprintf('E:/career/junior2/code/yc/dipolar/alpha%dN%dhx%dhz%ddtau02.mat',alpha,N_atom, hx, hz);
     
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
    Ut.matrix = { eye(2^N_atom) };
    U1 = expm( -i* Hz.matrix{1}* dtau );
    U2 = expm( -i* Hx.matrix{1}* dtau );
    
    loaded = false
end % if exist(filename,'file')

%% test operator 
Z0 = 1/N_atom * OperatorClass(N_atom,'z',1/2).matrix{1};
X0 = 1/N_atom * OperatorClass(N_atom,'x',1/2).matrix{1};

%% begin of for
ZZ_list = zeros(1, length(tau_list));
XX_list = zeros(1, length(tau_list));

for tau_index = 1:length(tau_list)

tau = tau_list(tau_index);
if loaded
    Ut = Ut_all{tau_index};
else
    Ut.matrix{1} = U1 * Ut.matrix{1} * U2;
    %Ut_all{end+1} = copy( Ut);
end

%% compute M(t -> \infty)
    Ut_t = Ut.matrix{1}^F_N ;
    ZZ_sum = 0;
    XX_sum = 0;
    for t = 1:F_N
        Ut_t = Ut.matrix{1} * Ut_t;
        ZZ_sum = ZZ_sum + trace(Ut_t' * Z0 * Ut_t * Z0);
        XX_sum = XX_sum + trace(Ut_t' * X0 * Ut_t * X0);
    end % for t = 1:IPR_N

    ZZ_list(tau_index) = real( ZZ_sum) / F_N/ 2^N_atom * 4*N_atom;
    XX_list(tau_index) = real( XX_sum) / F_N/ 2^N_atom * 4*N_atom;
    
end % for len_index = 1:length(N_atom_list)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end

%% plot log_IPR vs N_atom
ycPlot(24, tau_list, [ZZ_list; XX_list], '$J_z\tau$', '$<S(\infty)S(0)>_{\beta=0}$', ...,
    ['dipolar,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz), ...
    ' ,N_{atom}=',mat2str(N_atom), ' ,N_{period}=',mat2str(F_N),'$'], {'ZZ','XX'});

profile viewer