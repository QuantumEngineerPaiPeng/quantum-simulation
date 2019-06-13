%20190403
% compute fidelity of DQS
% Fid(t) = |U^t - U_tar^t|
% Fid1 = tr(U^dagger U_tar) / 2^N_atom = tr(UU) / 2^N_atom
% FidX = tr( U^d X U U_tar^d X U_tar ) / tr(XX), tr(XX) = N/4 * 2^N
% tr( U^d X U U_tar^d X U_tar ) = tr( UU X UU' X ), UU = U_tar U'
% of trotterized Ising chain with algebraically decaying
% interacions, also can be used for dipolar
% U = e^{-iHz\tau}e^{-iHx\tau}, U_tar=e^{-iH\tau}
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 4/3

%% runtime parameters
N_all = [8];

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;
Hx_len = 0;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = ['kp'];

dtau = 0.05;
tau_list = [dtau:dtau:3];
UU_tau = 1.;

Jt_list = [50];
n_list = Jt_list;

profile on
profile clear


%% begin of for diff N
[Fid, FidU, FidX, FidZ, XX_list] = deal( zeros(length(N_all), length(tau_list), length(n_list)) );
%Fid_pure = zeros(length(N_all), length(tau_list), length(n_list), psi_N);
for N_index = 1: length(N_all)
    
N_atom = N_all(N_index)
Jza = 1;% Jz * (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;
filename=sprintf('E:/career/junior2/code/yc/algebraic/alpha%dN%dhx%dhz%ddtau02.mat',alpha,N_atom, hx, hz);

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);
Hz = (1-Hx_len)* Hz;

Hx=hx* OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);
Hx = Hx + Hx_len* Hz;

H_tar = Hz + Hx;

Ut=OperatorClass(N_atom);
Ut.matrix = { speye(2^N_atom) };
Ut.symmetrize(sym);
U1 = H2U(Hz, dtau);
U2 = H2U(Hx, dtau);

%% test operator
X0=OperatorClass(N_atom,'x',1/2);
X0.symmetrize(sym);
Z0=OperatorClass(N_atom,'z',1/2);
Z0.symmetrize(sym);
Y0=OperatorClass(N_atom,'y',1/2);
Y0.symmetrize(sym);

%% test initial pure states
%{
y_pi_2 = H2U(OperatorClass(N_atom,'y',1/2), pi/2);
psi0 = y_pi_2.matrix{1} * initial_psi( 'up/down', N_atom, psi_N );
rho0 = cell(1,psi_N) ;
for psi0_index = 1: psi_N
    rho0{psi0_index} = OperatorClass(N_atom);
    rho0{psi0_index}.matrix = { psi0(:, psi0_index) * psi0(:, psi0_index)' };
    rho0{psi0_index}.symmetrize(sym);
end
%}

%% begin of for diff tau
for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    %dn = floor( Jt_end/ tau );
    
    Ut = U1 * Ut * U2;
    U_tar = H2U(H_tar, tau);
    
    if tau >= UU_tau && tau-dtau < UU_tau
        plot_UU = U_tar_t * U_t';
        1
    end

    for n_index=1: length(Jt_list)
        U_tar_t = U_tar^floor(Jt_list(n_index)/ tau);
        U_t = Ut^floor(Jt_list(n_index)/ tau);
        UU = U_tar_t * U_t';
        FidU(N_index, tau_index, n_index) = abs( trace(UU) )/ 2^N_atom;
        FidX(N_index, tau_index, n_index) = trace( UU* X0* UU'* X0 )/ 2^N_atom *4/ N_atom;
        FidZ(N_index, tau_index, n_index) = trace( UU* Z0* UU'* Z0 )/ 2^N_atom *4/ N_atom;
        FidY(N_index, tau_index, n_index) = trace( UU* Y0* UU'* Y0 )/ 2^N_atom *4/ N_atom;
        XX_list(N_index, tau_index, n_index) = abs( trace( U_t* X0* U_t'* X0 ) )/ 2^N_atom *4/ N_atom;
        %{ 
        decompose UU
        for psi0_index = 1: psi_N
            Fid_pure(N_index, tau_index, n_index, psi0_index) = abs( trace( UU* rho0{psi0_index} ) );
        end
        %}
        %{
        dU = U_t - U_tar_t;
        max_list = [];
        for kp_index = 1: length(dU.matrix)
            if size( dU.matrix{kp_index}, 1 )== 0
                continue
            end
            max_list(end+1) = max( abs( dU.matrix{kp_index} ) ,[],'all');
        end
        Fid(N_index, tau_index, n_index) = max( max_list);
        %}
    end

end % for tau_index = 1:length(tau_list)

end % for N_index = 1: length(N_all)

%% fix N. plot Fid(tau, t)
%{
extra = {{'x_par', pi*2/3/dtau}};
figure(1)
clf;

subplot(1,2,1)
ycPlot_matrix( ...,
    squeeze( FidU(1,:,:) )', '$ h_x\tau$', '$ n $' ...
    , ['Ising U fidelity,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),'$'], extra...,
)
x_n = length(tau_list);
xticks(linspace(1, x_n+1, 4) )
xticklabels({0,1,2,3})

subplot(1,2,2)
ycPlot_matrix( ...,
    squeeze( FidX(1,:,:) )', '$ h_x\tau$', '$ n $' ...
    , ['Ising X fidelity,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),'$'], extra...,
)
xticks(linspace(1, x_n+1, 4) )
xticklabels({0,1,2,3})
%}

%% fix n, N. plot Fid vs tau
n0_index = 1;
legends = {'U', 'X', 'Y', 'Z'};
extra = {{'x_par',  1.1}};

figure(3)
clf;

subplot(1,2,1)
ycPlot(tau_list, [squeeze( FidU(end, :, n0_index )); squeeze( FidX(end, :, n0_index)); ...
    squeeze( FidY(end, :, n0_index)); squeeze( FidZ(end, :, n0_index))], '$ h_x\tau$' ...
    , 'Fidelity', ['dipolar, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,Jt=',mat2str(Jt_list(n0_index)),'$'], legends, extra);
ylim([-0.2 1])
n0_index = 2;
subplot(1,2,2)
ycPlot(tau_list, [squeeze( FidU(end, :, n0_index )); squeeze( FidX(end, :, n0_index)); ...
    squeeze( FidY(end, :, n0_index)); squeeze( FidZ(end, :, n0_index))], '$ h_x\tau$' ...
    , 'Fidelity', ['dipolar, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,Jt=',mat2str(Jt_list(n0_index)),'$'], legends, extra);
ylim([-0.2 1])
%% relation between fidelity and XX correlation
figure(11)
clf;
legends = {'$tr(X(t)X_{id}(t))$', '$tr(X(t)X(0))$', '$tr(X(t)X(0))tr(X_{id}(t)X(0))$'};
ycPlot(tau_list, [squeeze( FidX(end, :, end));XX_list(end, :, end); XX_list(end, :, end)*XX_list(end, 1, end)] ...
    , '$ h_x\tau$', 'value', ['dipo, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,Jt=',sprintf('%1.0e',Jt_list(end)),'$'], legends, extra);
%ylim([-0.2 1])

%% fix n. plot Fid vs tau for diff N
n0_index = 1;% find(n_list>=n0, 1);
legends = cell(0);
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end
extra = {{'x_par',  1.2}};

figure(88)
clf;
ycPlot(tau_list, [squeeze( FidX(:, :, n0_index ))], '$ h_x\tau$' ...
    , '$F_X$', ['dipo,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    ,' ,Jt=',sprintf('%1.0e',Jt_list(end)),'$'], legends, extra);

%{
%% fix n,N. plot FidX vs tau for diff psi0
n0 = 50;
n0_index = find(n_list>=n0, 1);
legends = cell(0);
N_atom = N_all(end);
for psi0_index = 1:psi_N
    legends{end+1} = ['$<\psi_0| 2X |\psi_0>=',mat2str( N_atom+1-psi0_index ),'$'];
end
extra = {{'x_par',  pi*2/3}};

figure(4)
clf;
ycPlot(tau_list, [squeeze( Fid_pure(end, :, n0_index,: ))'], '$ h_x\tau$' ...
    , '$F_X$', ['Ising,$\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom), ' ,n_{period}=',mat2str(n0),'$'], legends, extra);


%% fix N, n, tau, plot matrix UU
%extra = {{'x_par', pi*2/3/dtau}};
figure(5)
clf;
plot_U = plot_UU.matrix{1};
subplot(1,2,1)
ycPlot_matrix( ...,
    abs(plot_U), 'z-spin basis', 'z-spin basis' ...
    , ['Ising, $|{U^n}^\dagger U^n_0|, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], {}...,
)
x_n = 2^N_atom;
xticks(linspace(1, x_n+1, 2) )
xticklabels({'$|\uparrow\cdots \uparrow\rangle$','$|\downarrow\cdots \downarrow\rangle$'})
yticks(linspace(1, x_n+1, 2) )
yticklabels({'$|\uparrow\cdots \uparrow\rangle$','$|\downarrow\cdots \downarrow\rangle$'})


subplot(1,2,2)
plot_U = (y_pi_2 * plot_UU * y_pi_2');
sort_key = [];
for psi_index = 0: 2^N_atom -1
    sort_key(end+1) = length(find( dec2bin(psi_index)=='1'));
end
plot_U = sortrows( [sort_key; plot_U.matrix{1}]', 1);
plot_U = plot_U(:, 2:end)';
plot_U = sortrows( [sort_key; plot_U']', 1);
plot_U = plot_U(:, 2:end);

ycPlot_matrix( ...,
    abs(plot_U), '$ 2X$', '$ 2X$' ...
    , ['$N_{atom}=',mat2str(N_atom), ' ,h_x\tau=',mat2str(UU_tau),' ,n=',mat2str(n0),'$'], {}...,
)
%caxis([-log(2^N_atom) 0]);
x_n = 2^N_atom;
xticks(linspace(1, x_n+1, 5) )
xticklabels({mat2str(-N_atom), mat2str(-N_atom/4), 0, mat2str(N_atom/4), mat2str(N_atom)})
yticks(linspace(1, x_n+1, 5) )
yticklabels({mat2str(-N_atom), mat2str(-N_atom/4), 0, mat2str(N_atom/4), mat2str(N_atom)})
%}

profile viewer