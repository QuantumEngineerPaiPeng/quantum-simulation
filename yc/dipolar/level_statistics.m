%20190520
% level statistics of
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/20

%% runtime parameters
N_all = [6:2:8];

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;
hz_dis = 1;

model = [1, -1/2, -1/2];
bc = 'p';
sym = [];

dtau = 0.05;
tau_list = [0:dtau:3];

F_N = [10^2];
repeat = 10;

load_save = false;
diago = true;
    
profile on
profile clear

%% begin of for diff N
[r_list] = deal( zeros(repeat, length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
for rep_index = 1: repeat
    
N_atom = N_all(N_index)
Jza = Jz; %* (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%%  construct Hamiltonian and Floquet
disorder_fields = (2*rand(1, N_atom)-1) * hz_dis;
disorder_fields = [disorder_fields; zeros(2, N_atom)];
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'x',1/2) ...
    + 1/2*randOpe(N_atom, disorder_fields) ; % if disorder, should all be 'z'
Hz.symmetrize(sym);
Hx=hx* OperatorClass(N_atom,'z',1/2);
Hx.symmetrize(sym);

% find the largest subspace
%{
blockN = length(Hz.matrix);
max_index = 0;
max_dim = 0;
for block_index = 1: blockN
    if size(Hz.matrix{block_index},1) > max_dim
        max_dim = size(Hz.matrix{block_index}, 1);
        max_index = block_index;
    end
end
%}
if hz_dis > 0.01
    max_index = 1;
else
    max_index = 2;
end

if load_save && exist(filename,'file')
    load(filename);
    loaded = true;
else
    Ut_all = cell(0);
    Ut = speye(size(Hz.matrix{max_index} ,2));
    U1 = expm(-i*Hz.matrix{max_index}*dtau);
    U2 = expm(-i*Hx.matrix{max_index}*dtau);
    loaded = false;
end % if exist(filename,'file')

%% begin of for diff tau
for tau_index = 1:length(tau_list)

tau = tau_list(tau_index);
if tau_index == 0
    pse_e = sort( angle( eig(Hz.matrix{max_index}+Hx.matrix{max_index}) ) );
    s_list = pse_e(2:end) - pse_e(1:end-1);
    this_r_list = [];
    for eig_index = 1: length(s_list)-1
        this_r_list(end+1) = min(s_list(eig_index), s_list(eig_index+1))/ max(s_list(eig_index), s_list(eig_index+1));
    end
    r_list(rep_index, N_index, tau_index) = mean(this_r_list);
    continue
end

Ut = U1 * Ut * U2;

%Ut.diagonalize();
pse_e = sort( angle( eig(Ut) ) );
s_list = pse_e(2:end) - pse_e(1:end-1);
this_r_list = [];
for eig_index = 1: length(s_list)-1
    this_r_list(end+1) = min(s_list(eig_index), s_list(eig_index+1))/ max(s_list(eig_index), s_list(eig_index+1));
end
r_list(rep_index, N_index, tau_index) = mean(this_r_list);

end % for tau_index = 1:length(tau_list)

end % for rep_index = 1: repeat
end % for N_index = 1: length(N_all)

%% finite size scaling
legends = cell(0);
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end
%legends{end+1} = 'analytic';

figure(26)
clf;
%{
subplot(1,2,1);
ycPlot(tau_list, ZZ_list, '$J_z\tau$', '$<H_{dipz}(\infty)H_{dipz}(0)>_{\beta=0}$', ...,
    ['dipolar, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{period}=\infty$'], legends);
subplot(1,2,2);
%}
mean( r_list(:,1:end,:), 1)
sqrt( var( r_list(:,1:end,:), 1) )
extra = {{'y_par',0.386}, {'y_par',0.536}};
ycPlot(tau_list, squeeze( mean( r_list(:,1:end,:), 1) ), '$h_x\tau$', '$\langle r\rangle$', ...,
    ['dipo, symmetric sector, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz), '$'], legends, extra);
%ylim([0.48 0.54])
profile viewer