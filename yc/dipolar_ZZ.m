function [ZZ_list, XX_list] = dipolar_ZZ( ...,
    N_atom, hx, Jz, hz ...
    ,dtau, tau_list, t_intervals, extra ...
)
%20190319 
% also can compute dipo(t)dipo(0)
% compute observable ZZ(t) and XX(t) of dipolar_z with hx and hz
% interacions at T = \infty, sometimes also ising
% ZZ_list(length(t_intervals), length(tau_list))
% t_intervals = { [t11, t12], [t21, t22], ..., 'ED' }
%
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N *4/N, Z_0 = S^z
% normalize by ZZ(0) = 1/(4N)
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 4/28, add disorder

%% fixed parameters
load_save = false;
ising = false;
disorder = false;
alpha = 100;
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'load_save')
        load_save = extra_item{2};
    end
    if strcmp( extra_item{1}, 'ising')
        ising = extra_item{2};
    end
    if strcmp( extra_item{1}, 'alpha')
        alpha = extra_item{2};
    end
    if strcmp( extra_item{1}, 'hz_dis')
        disorder = true;
        hz_dis = extra_item{2};
    end
end
if disorder == false
    hz_dis = 0;
end

%% model
if ising == false
    %alpha = 3;
    model = [1, -1/2, -1/2];
else
    %alpha = 100;
    model = [1, 0, 0];
end
Jza = (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) ); %Jz;
cplist = Jza * [1:N_atom-1].^(-alpha) ;

bc = 'p';
if ~disorder
    sym = ['kpz'];
else
    sym = [];
end

filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau02.mat',N_atom, hx, hz);

%%  construct Hamiltonian and Floquet
disorder_fields = (2*rand(1, N_atom)-1) * hz_dis;
disorder_fields = [zeros(2, N_atom); disorder_fields];
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'x',1/2) ...
    + 1/2*randOpe(N_atom, disorder_fields) ; % if disorder, should all be 'z'
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'z',1/2);
Hx.symmetrize(sym);

if load_save && exist(filename,'file')
    load(filename);
    loaded = true;
else
    Ut_all = cell(0);
    Ut=OperatorClass(N_atom);
    Ut.matrix = { speye(2^N_atom) };
    Ut.symmetrize(sym);
    U1 = H2U(Hz, dtau);
    % U1 = H2U(Hz+2*H_dis, dtau/2) * H2U(Hz-H_dis, dtau/2);
    U2 = H2U(Hx, dtau);
    
    loaded = false;
end % if exist(filename,'file')

%% test operator 
Z0 = OperatorClass(N_atom,model,1/4,bc,cplist);
X0 = OperatorClass(N_atom,'z',1/2);
Z0.symmetrize(sym);
X0.symmetrize(sym);

%% begin of for
ZZ_list = zeros(length(t_intervals), length(tau_list));
XX_list = zeros(length(t_intervals), length(tau_list));

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    if loaded
        Ut = Ut_all{tau_index};
    else
        Ut = U1 * Ut * U2;
        %Ut_all{end+1} = copy( Ut);
    end
    
    for t_index = 1: length(t_intervals)
        temp =  ZZ_longtime(Ut, floor( t_intervals{t_index}/ tau ), {Z0,X0});
        ZZ_list(t_index, tau_index) = temp(1, :)*4/N_atom;
        XX_list(t_index, tau_index) = temp(2, :)*4/N_atom;
    
    end % for t_index = 1: length(t_list)

end % for tau_index = 1:length(tau_list)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end


end





