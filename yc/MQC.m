function [I_qtz, I_qtx] = MQC( ...,
    N_atom, hx, Jz, hz ...
    ,dtau, tau_list, t_intervals, load_save, ising ...
)
%20190323 
% compute observable I_q(t) of dipolar_z with hx and hz
% interacions at T = \infty, sometimes also ising
% I_qt(N_atom+1 ,length(t_intervals)), I_qt{t_index, tau_Index}(q, t)
% t_intervals = { [t11, t12], [t21, t22], ..., 'ED' }
%
% S_mt = S_m(t) = 2^(-L) tr( O0 U(m,t)' O0 U(m,t) )
% U(m, t) = U_t^(-t) exp(im pi PZ/L) U_t^t
% H = Hz + Hx
% Hz = Jz \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/23

%% fixed parameters
if nargin == 8 | ising == false
    ising = false;
    alpha = 3;
    model = [-1/2, -1/2, 1];
elseif ising == true
    alpha = 100;
    model = [0, 0, 1];
end
Jza = Jz;
cplist = Jza * [1:N_atom-1].^(-alpha) ;
m_list = [0: 2*N_atom-1];

bc = 'p';
% sym = ['kp'];
sym = [];

filename=sprintf('E:/career/junior2/code/yc/dipolar/N%dhx%dhz%ddtau02.mat',N_atom, hx, hz);

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);

Hx=hx* OperatorClass(N_atom,'x',1/2);
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
    U2 = H2U(Hx, dtau);
    
    loaded = false;
end % if exist(filename,'file')

%% test operator 
Z0 = OperatorClass(N_atom,'z',1/2);
X0 = OperatorClass(N_atom,'x',1/2);
Z0.symmetrize(sym);
X0.symmetrize(sym);

%% P operator
PZ = Z0;

%% begin of for
[Z_MQC, X_MQC, I_qtz, I_qtx] = deal( cell(length(t_intervals), length(tau_list)) ) ;

for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    if loaded
        Ut = Ut_all{tau_index};
    else
        Ut = U1 * Ut * U2;
        Ut_all{end+1} = copy( Ut);
    end
    
    for t_index = 1: length(t_intervals)
        t_intervals{t_index}(1)/ tau_list(tau_index)
        t_list = floor( [t_intervals{t_index}(1)/ tau_list(tau_index): ...
            t_intervals{t_index}(2)/ tau_list(tau_index)] );
        if length(t_list) > 100
            t_list = t_list(1:1+floor(length(t_list)/100) :end);
        end
        t_list
        Z_MQC{t_index,tau_index} =  MQC_S(Ut_all{tau_index}, Z0, m_list, t_list, Z0);
        X_MQC{t_index,tau_index} =  MQC_S(Ut_all{tau_index}, X0, m_list, t_list, Z0);
    
    end % for t_index = 1: length(t_list)

end % for tau_index = 1:length(tau_list)
assignin('base','Z_M',Z_MQC)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end

%% do fourier transform to get MQC intensity
for t_index = 1: length(t_intervals)
    for tau_index = 1: length(tau_list)
        I_qtz{t_index, tau_index} = real( ifft(Z_MQC{t_index, tau_index}, [], 1) );
        I_qtx{t_index, tau_index} = real( ifft(X_MQC{t_index, tau_index}, [], 1) );
    end
end

end





