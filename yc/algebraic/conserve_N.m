%20190418
% finite size scaling of:
% investigate the approximate conservation of X in traversed Ising chain
% or Ising chain with algebraically decaying
% interacions, also can be used for dipolar
% also for such models under driven, one driven term: H_1 = h_x X 
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 4/18

%% runtime parameters
N_list = [8:12];

alpha = 100;
hx = 1;
Jz = 1;
hz = 0;

tau = 2;
Hx_len = 0.;

model = [-1/2,-1/2, 1];
bc = 'p';
sym = [];

n_start = 10^4;
n_ave = 50;
whole_U = false;
    
profile on
profile clear

%% begin of diff size
%Un = cell(1, length(n_list));
U_diag_list = cell(2, length(N_list));
U_diag_err = cell(1, length(N_list));
X_loc = cell(1, length(N_list));

for N_index = 1: length(N_list)
    N_atom = N_list(N_index)

Jza = Jz;
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);
Hz = (1-Hx_len)* Hz;

Hx=hx* OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);
Hx = Hx + Hx_len* Hz;

if abs(tau) > 1e-10
    Un1 = H2U(Hx, tau)* H2U(Hz, tau);
else
    H_tar = Hz + Hx;
    Un1 = H2U(H_tar, 1);
end

sumU_list = zeros(N_atom+1, N_atom+1, n_ave);
for n_index=1: n_ave
    if n_index == 1
        Un = Un1^n_start;
    else
        Un = Un * Un1;
    end
    sumU_list(:,:, n_index) = sumU_X(Un);
end
sumU = mean( sumU_list, 3);
sumU_err = sqrt(var( sumU_list, 0, 3 ));
sumU_n1 = sumU_X(Un1);

U_diag_list{1, N_index} = diag( sumU_n1 );
U_diag_err{ N_index} = diag( sumU_err );
U_diag_list{2, N_index} = diag( sumU );
max(diag( sumU_err ))

for nd = 0:N_atom
    temp = [0 0];
    for nd_c = 0:N_atom
        temp(1) = temp(1) + (nd-nd_c)^2 * sumU_n1(nd+1, nd_c+1);
        temp(2) = temp(2) + (nd-nd_c)^2 * sumU(nd+1, nd_c+1);
    end
    X_loc{1, N_index}(nd+1) = sqrt(temp(1))/ N_atom;
    X_loc{2, N_index}(nd+1) = sqrt(temp(2))/ N_atom;
end

end % for N_index
%% plot sum of abs( U_tar^n) in X basis
if model(1)==0
    model_ = 'Ising'
else
    model_ = 'dipolar'
end

figure(16)
clf;
legends = cell(0);
for N_index = 1:length(N_list)
    legends{end+1} = ['$N_{atom}=',mat2str( N_list(N_index) ),'$'];
end

subplot(1,2,1)
hold on
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), U_diag_list{1,N_index}, '$2X/N_{atom}$', '$\sum_{jk}|P_X U^n P_X|_{jk}^2$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=1$'], legends);
end

subplot(1,2,2)
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), U_diag_list{2,N_index}, '$2X/N_{atom}$', '$\sum_{jk}|P_X U^n P_X|_{jk}^2$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=',mat2str(n_start),'$'], legends);
end

%% plot X_loc/N
figure(17)
clf;
legends = cell(0);
for N_index = 1:length(N_list)
    legends{end+1} = ['$N_{atom}=',mat2str( N_list(N_index) ),'$'];
end

subplot(1,2,1)
hold on
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), X_loc{1,N_index}, '$2X/N_{atom}$', '$X_{loc}/N_{atom}$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=1$'], legends);
end

subplot(1,2,2)
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), X_loc{2,N_index}, '$2X/N_{atom}$', '$X_{loc}/N_{atom}$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=', mat2str(n_start),'$'], legends);
end

profile viewer


%% plot X_loc
figure(18)
clf;
legends = cell(0);
for N_index = 1:length(N_list)
    legends{end+1} = ['$N_{atom}=',mat2str( N_list(N_index) ),'$'];
end

subplot(1,2,1)
hold on
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), X_loc{1,N_index}* N_list(N_index), '$2X/N_{atom}$', '$X_{loc}$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=1$'], legends);
end

subplot(1,2,2)
for N_index = 1: length(N_list)
    ycPlot(linspace(-1,1, N_list(N_index)+1), X_loc{2,N_index}* N_list(N_index), '$2X/N_{atom}$', '$X_{loc}$', ...,
        [model_,', $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
        , ' ,\tau=',mat2str(tau), ' ,n=', mat2str(n_start),'$'], legends);
end

profile viewer

%% function of \sum_{jk}|P_X U_0^n P_X|_{jk}^2
function sumU = sumU_X(Un)

N_atom = Un.L;
%% change to X basis
y_pi_2 = H2U(OperatorClass(N_atom,'y',1/2), pi/2);
sort_key = [];
for psi_index = 0: 2^N_atom -1
    sort_key(end+1) = length(find( dec2bin(psi_index)=='1'));
end
    temp_U = (y_pi_2 * Un * y_pi_2');
    
    temp_U = sortrows( [sort_key; temp_U.matrix{1}]', 1);
    temp_U = temp_U(:, 2:end)';
    temp_U = sortrows( [sort_key; temp_U']', 1);
    UnX = temp_U(:, 2:end);

%% sum for diff X block
key_sorted = sort(sort_key);
X_begin_list = zeros(1, N_atom+2);
for n_down = 0: N_atom
    X_begin_list(n_down+1) = find(key_sorted>=[n_down], 1);
end
X_begin_list(end) = 2^N_atom+1;

sumU = zeros(N_atom+1, N_atom+1);
    for nd_r = 0:N_atom
        for nd_c = 0:N_atom
            sumU(end-nd_r, end-nd_c) = sum(sum( abs(UnX( ...
                X_begin_list(nd_r+1):X_begin_list(nd_r+2)-1, ...
                X_begin_list(nd_c+1):X_begin_list(nd_c+2)-1 )  ).^2 )) ./ ...
                sqrt( (X_begin_list(nd_r+2)-X_begin_list(nd_r+1))* ...
                (X_begin_list(nd_c+2)-X_begin_list(nd_c+1)));
        end
    end % for nd_r = N_atom: -2: 0

end

%% function of matrix plot
function arg_out = approx_conserve_plot(ma_in, ma_in_X, N_atom, hx, hz, plot_n)
    
subplot(1,2,1)
ycPlot_matrix( ...,
    abs(ma_in.matrix{1}), 'z-spin basis', 'z-spin basis' ...
    , ['Ising, $|U^n_0|, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], {}...,
)
xticks(linspace(1, x_n+1, 2) )
xticklabels({'$|\uparrow\cdots \uparrow\rangle$','$|\downarrow\cdots \downarrow\rangle$'})
yticks(linspace(1, x_n+1, 2) )
yticklabels({'$|\uparrow\cdots \uparrow\rangle$','$|\downarrow\cdots \downarrow\rangle$'})


subplot(1,2,2)
ycPlot_matrix( ...,
    abs(ma_in_X), '$ 2X$', '$ 2X$' ...
    , ['$N_{atom}=',mat2str(N_atom),' ,n=',mat2str(plot_n),'$'], {}...,
)
%caxis([-log(2^N_atom) 0]);
x_n = 2^N_atom;
xticks(linspace(1, x_n+1, 5) )
xticklabels({mat2str(-N_atom), mat2str(-N_atom/4), 0, mat2str(N_atom/4), mat2str(N_atom)})
yticks(linspace(1, x_n+1, 5) )
yticklabels({mat2str(-N_atom), mat2str(-N_atom/4), 0, mat2str(N_atom/4), mat2str(N_atom)})
end





