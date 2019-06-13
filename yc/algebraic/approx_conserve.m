%20190406
% investigate the approximate conservation of X in traversed Ising chain
% or Ising chain with algebraically decaying
% interacions, also can be used for dipolar
% also for such models under driven, one driven term: H_1 = h_x X 
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 4/6

%% runtime parameters
N_atom = 10;

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;

trotter = true;
tau = 1;
Hx_len = 0.;

model = [-1/2, -1/2, 1];
bc = 'p';
sym = [];

n_list = [ 10.^[0:2], 10^4+100 ];
Un = cell(1, length(n_list));
whole_U = false;
    
profile on
profile clear

Jza = Jz;
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%%  construct Hamiltonian and Floquet
Hz=OperatorClass(N_atom,model,1/4,bc,cplist) + hz* OperatorClass(N_atom,'z',1/2);
Hz.symmetrize(sym);
Hz = (1-Hx_len)* Hz;

Hx=hx* OperatorClass(N_atom,'y',1/2);
Hx.symmetrize(sym);
Hx = Hx + Hx_len* Hz;

if trotter
    U0 = H2U(Hx, tau)* H2U(Hz, tau) ;
    for n_index=1: length(n_list)
        Un{n_index} = U0^n_list(n_index);
    end
else
    tau = 0;
    H_tar = Hz + Hx;

    for n_index=1: length(n_list)
        Un{n_index} = H2U(H_tar, n_list(n_index));
    end
end

%% change to X basis
y_pi_2 = H2U(OperatorClass(N_atom,'y',1/2), pi/2);
sort_key = [];
for psi_index = 0: 2^N_atom -1
    sort_key(end+1) = length(find( dec2bin(psi_index)=='1'));
end
UnX = cell(1, length(Un));
for n_index = 1: length(n_list)
    temp_U = (y_pi_2 * Un{n_index} * y_pi_2');
    
    temp_U = sortrows( [sort_key; temp_U.matrix{1}]', 1);
    temp_U = temp_U(:, 2:end)';
    temp_U = sortrows( [sort_key; temp_U']', 1);
    UnX{n_index} = temp_U(:, 2:end);
end

%% sum for diff X block
key_sorted = sort(sort_key);
X_begin_list = zeros(1, N_atom+2);
for n_down = 0: N_atom
    X_begin_list(n_down+1) = find(key_sorted>=[n_down], 1);
end
X_begin_list(end) = 2^N_atom+1;

sumU = zeros(N_atom+1, N_atom+1, length(n_list));
for n_index = 1: length(n_list)
    for nd_r = 0:N_atom
        for nd_c = 0:N_atom
            sumU(end-nd_r, end-nd_c, n_index) = sum(sum( abs(UnX{n_index}( ...
                X_begin_list(nd_r+1):X_begin_list(nd_r+2)-1, ...
                X_begin_list(nd_c+1):X_begin_list(nd_c+2)-1 )  ).^2 )) ./ ...
                sqrt( (X_begin_list(nd_r+2)-X_begin_list(nd_r+1))* ...
                (X_begin_list(nd_c+2)-X_begin_list(nd_c+1)));
        end
    end % for nd_r = N_atom: -2: 0
end

sumU_diag = zeros(N_atom+1, length(n_list));
for n_index = 1: length(n_list)
    for nd_r = 0:N_atom
        sumU_diag(nd_r+1, n_index) = sumU(nd_r+1, nd_r+1, n_index);
    end % for nd_r = N_atom: -2: 0
end

%% plot matrix U_tar^n, in Z and X basis
if whole_U == true
for n_index = 1: length(n_list)
    figure(10 + n_index)
    clf;
    approx_conserve_plot(Un{n_index},UnX{n_index}, N_atom, hx, hz, n_list(n_index));
end
end
%% plot sum of abs( U_tar^n) in X basis

figure(15)
clf;
subplot(1,2,1)
n_index = 4;
ycPlot_matrix( ...,
    sumU(:,:,n_index), '$X_1$', '$X_2$' ...
    , ['$\sum_{jk}|P_{X_1} U_0^n P_{X_2}|_{jk}^2/\sqrt{d(X_1)\cdot d(X_2)}, n=', mat2str(n_list( n_index) )...
    , '$'], {}...,
)
xticks(linspace(1,N_atom+1, 3) )
xticklabels({mat2str(-N_atom/2), 0, mat2str(N_atom/2)})
yticks(linspace(1, N_atom+1, 3) )
yticklabels({mat2str(-N_atom/2), 0, mat2str(N_atom/2)})

legends = cell(0);
for n_index = 1:length(n_list)
    legends{end+1} = ['$n=',mat2str( n_list(n_index) ),'$'];
end
subplot(1,2,2)
ycPlot([-N_atom/2:N_atom/2], sumU_diag, '$X$', '$\sum_{jk}|P_X U_0^n P_X|_{jk}^2$', ...,
    ['dipolar, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,\tau=',mat2str(tau), ' ,N_{atom}=',mat2str(N_atom),'$'], legends);

profile viewer



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





