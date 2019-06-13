%20190530
% level statistics of
% H = Hz + Hx, averaged over all sectors
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/30

%% runtime parameters
N_all = [10:2:12];

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;
hz_dis = 0;

model = [1, -1/2, -1/2];
bc = 'p';
sym = ['kpz'];

dtau = 0.05;
tau_list = [0:dtau:3];

F_N = [10^2];
repeat = 1;
min_dim = 50;

load_save = false;
diago = true;
    
profile on
profile clear

%% begin of for diff N
[r_list, dr_list] = deal( zeros(repeat, length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
for rep_index = 1: repeat
    
N_atom = N_all(N_index)
Jza = Jz; %* (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%%  construct Hamiltonian and Floquet
disorder_fields = (2*rand(3, N_atom)-1) * hz_dis;
%disorder_fields = [zeros(2, N_atom); disorder_fields];
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'x',1/2) ...
    + 1/2*randOpe(N_atom, disorder_fields) ; % if disorder, should all be 'z'
Hz.symmetrize(sym);
Hx=hx* OperatorClass(N_atom,'z',1/2);
Hx.symmetrize(sym);

% find all sectors that are large enough
blockN = length(Hz.matrix);
sector_list = [];
for block_index = 1: blockN
    if size(Hz.matrix{block_index},1) > min_dim
        sector_list(end+1) = block_index;
    end
end
sector_list

Hz.diagonalize();
Hx.diagonalize();
U1 = H2U(Hz, dtau);
U2 = H2U(Hx, dtau);

%% begin of for diff tau
for tau_index = 1:length(tau_list)
    tau = tau_list(tau_index);
    if tau_index == 1
        Ut = Hz+ Hx;
    elseif tau_index == 2
        Ut = U1* U2;
    else
        Ut = U1 * Ut * U2;
    end

    Ut.diagonalize();
    
    mean_r_list = [];
    for sector_index = 1: length( sector_list)
        this_r_list = [];
        pse_e = sort( angle( Ut.eigsys{sector_list(sector_index)}.D ) );
        s_list = pse_e(2:end) - pse_e(1:end-1);
        for eig_index = 1: length(s_list)-1
            this_r_list(end+1) = min(s_list(eig_index), s_list(eig_index+1))/ max(s_list(eig_index), s_list(eig_index+1));
        end
        mean_r_list(end+1) = mean( this_r_list);
    end % for sector_index = 1: sector_list
    r_list(rep_index, N_index, tau_index) = mean( mean_r_list);
    dr_list(rep_index, N_index, tau_index) = sqrt( var( mean_r_list) / length(sector_list) );

end % for tau_index = 1:length(tau_list)

end % for rep_index = 1: repeat
end % for N_index = 1: length(N_all)

%% finite size scaling
legends = cell(0);
for N_index = 1:length(N_all)-1
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end
legends{end+1} = ['$N_{atom}=18$, 3 sectors'];
%legends{end+1} = ['$N_{atom}=18$, 10 sectors'];
%legends{end+1} = ['$N_{atom}=20$, 1 sector'];

figure(26)
clf;
%{
subplot(1,2,1);
ycPlot(tau_list, ZZ_list, '$J_z\tau$', '$<H_{dipz}(\infty)H_{dipz}(0)>_{\beta=0}$', ...,
    ['dipolar, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{period}=\infty$'], legends);
subplot(1,2,2);
%}
extra = {{'y_par',0.386}, {'y_par',0.536} ,{'err', squeeze( var( r_list, 1) )}};
ycPlot(tau_list, squeeze( mean( r_list(:,:,:) ,1) ), '$h_x\tau$', '$\langle r\rangle$', ...,
    ['Ising, 10 sector ave, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz), '$'], legends, extra);
ylim([0.42 0.54])
%yticks([0.38:0.02:0.54])
profile viewer