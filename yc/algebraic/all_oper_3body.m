%20190605
% correlations among opers: Cij = tr( opers{i}(t) opers{j} ), then
% diagonaize C
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 6/5

%% runtime parameters
N_all = [6:2:10];

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;

model = [1, -1/2, -1/2];
bc = 'p';
sym = [];

dtau = 0.1;
tau_list = [dtau:dtau:3];

oper_N = 1;
t_interval = 'ED'; % only 2 element
    
profile on
profile clear

%% begin of for diff N
OO_list = deal( zeros( oper_N,oper_N, length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
N_atom = N_all(N_index)
Jza = Jz; %* (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%% test operator 
X0 = OperatorClass(N_atom,'z',1/2);
X0.symmetrize(sym);
temp = 0;
opers = {}; 

sigma_list = [
    3,3,3; % XXX
    1,3,1; % ZXZ
    2,3,2; % YXY
    3,2,2; % XYY+YYX
    3,1,1; % XZZ+ZZX
    1,2,3; % ZYX+XYZ
    3,1,2; % XZY+YZX
    1,3,2; % ZXY+YXZ
];

for oper_index = 1: 3
    sig_list = sigma_list(oper_index,:);
    opers{end+1} = OperatorClass(N_atom);
    opers{end}.matrix = { sparse(2^N_atom, 2^N_atom) };
    for atom_index = 1: N_atom
        opers{end}.matrix{1} = opers{end}.matrix{1} + Kron3body(N_atom, sig_list(1), sig_list(2), sig_list(3) ...
            ,atom_index, mod(atom_index,N_atom)+1, mod(atom_index+1,N_atom)+1);
        if sig_list(1) ~= sig_list(3)
            opers{end}.matrix{1} = ( opers{end}.matrix{1} + Kron3body(N_atom, sig_list(3), sig_list(2), sig_list(1) ...
            ,atom_index, mod(atom_index,N_atom)+1, mod(atom_index+1,N_atom)+1) )/sqrt(2);
        end
    end
    opers{end}.matrix{1} = opers{end}.matrix{1} /8*4;
    opers{end}.symmetrize(sym);
end
opers = {(opers{2}+ opers{3})*sqrt(1/2)}

dipz = OperatorClass(N_atom,model,1/4,bc,cplist);
dipz.symmetrize(sym);

%% test operator: current J = XY-YX (when dipo is Z), here = YZ-ZY. energy current JE = \sum -i [hi_1, hi]
JP = OperatorClass(N_atom);
JP.matrix = { sparse(2^N_atom, 2^N_atom) };
for atom_index = 1: N_atom
    JP.matrix{1} = JP.matrix{1} + Kron2body(N_atom,2, 3,atom_index, mod(atom_index,N_atom)+1) ...
        - Kron2body(N_atom,3,2,atom_index, mod(atom_index,N_atom)+1);
end
JP.symmetrize(sym);

JE = OperatorClass(N_atom);
JE.matrix = { sparse(2^N_atom, 2^N_atom) };
for atom_index = 1: N_atom
    hi_1 = Kron2body(N_atom,1,1,atom_index, mod(atom_index,N_atom)+1) ...
        -1/2* Kron2body(N_atom,2,2,atom_index, mod(atom_index,N_atom)+1) ...
        -1/2* Kron2body(N_atom,3,3,atom_index, mod(atom_index,N_atom)+1);
    hi = Kron2body(N_atom,1,1,mod(atom_index,N_atom)+1, mod(mod(atom_index,N_atom)+1,N_atom)+1) ...
        -1/2* Kron2body(N_atom,2,2,mod(atom_index,N_atom)+1, mod(mod(atom_index,N_atom)+1,N_atom)+1) ...
        -1/2* Kron2body(N_atom,3,3,mod(atom_index,N_atom)+1, mod(mod(atom_index,N_atom)+1,N_atom)+1);
    JE.matrix{1} = JE.matrix{1} -i/4* (hi_1*hi -hi*hi_1);
end
JE.symmetrize(sym);

%%  construct Hamiltonian and Floquet

U1 = H2U(dipz, dtau);
U2 = H2U(hx* X0, dtau);

%% begin of for diff tau
for tau_index = 1:length(tau_list)

tau = tau_list(tau_index);
if tau < 0.001
    Ut = H2U( Hz+ Hx, 1);
elseif tau_index == 1
    '1st tau'
    Ut = U1* U2;
else
    Ut = U1 * Ut * U2;
end

extra = {{ 'cross', true }};
temp =  ZZ_longtime(Ut, t_interval, opers, extra);
OO_list(:,:, N_index,tau_index) = temp*4/N_atom;

end % for tau_index = 1:length(tau_list)

end % for N_index = 1: length(N_all)

%% diagonalize O_i O_j
OO_eig = zeros( oper_N, length(N_all), length(tau_list));
OO_v = OO_list;
for N_index = 1: length(N_all)
    for tau_index = 1: length(tau_list)
        [OO_v(:, :, N_index, tau_index), temp]= eig( OO_list(:,:, N_index, tau_index) );
        OO_eig(:, N_index, tau_index) = diag(temp);
    end
end
%save(['temp_mat.mat'],'N_all','tau_list','Jz','hx','hz','OO_list','OO_eig')

%% OO_v change with N_atom
%{
[N_dot, N_dot_ideal] = deal( zeros(4, 4, length(tau_list)) );
V_ideal = [0, 1/sqrt(3), -sqrt(2/3); -1/sqrt(2), 1/sqrt(3), 1/sqrt(6);  1/sqrt(2), 1/sqrt(3), 1/sqrt(6)];
for tau_index = 1: length(tau_list)
    N_dot_ideal(:, :,tau_index) = abs( abs( OO_v(:,:, end, tau_index)' * OO_v(:,:, end-1, tau_index) )- eye(oper_N) );
    %N_dot_ideal(:, :,tau_index) = abs( abs( OO_v(:,:, end, tau_index)' * V_ideal )- eye(oper_N) );
end

figure(30)
clf;

subplot(1,2,1)
ycPlot(tau_list, [squeeze( N_dot_ideal(1,1,:)), squeeze( N_dot_ideal(2,2,:)), squeeze( N_dot_ideal(3,3,:))...
    , squeeze( N_dot_ideal(4, 4,:))] ...
    , '$h_x\tau$', '$1- |v_i \cdot v_i^0|$', ...,
    ['eig vectors $v_i$ of $O_\mu O_\nu, v^0 = ZZ-YY, Hei, Dip_X','$'], {'1','2','3','4'});

%V_ideal = [1/sqrt(2), 1/sqrt(3), 1/sqrt(6); -1/sqrt(2), 1/sqrt(3), 1/sqrt(6);  0, 1/sqrt(3),-sqrt(2/3) ];
for tau_index = 1: length(tau_list)
    N_dot_ideal(:, :,tau_index) = abs( abs( OO_v(:,:, end, tau_index)' * V_ideal )- eye(oper_N) );
end
subplot(1,2,2)
ycPlot(tau_list, [squeeze( N_dot_ideal(1,1,:)), squeeze( N_dot_ideal(2,2,:)), squeeze( N_dot_ideal(3,3,:))] ...
    , '$h_x\tau$', '$1- |v_i \cdot v_i^0|$', ...,
    ['$Dip_z+X, \frac{h_x}{J_z}=',mat2str(hx) ,', v^0 = XX-YY, Hei, Dip_Z','$'], {'1','2','3'});
%}

%% finite size scaling
legends = cell(0);
%legends{end+1} = 'analytic';

figure(27)
clf;

subplot(1,2,1)
hold on
line_type_list = {'-','--',':', '-.'};
[colors, lines] = deal( {} );
legend_lines = [];
for N_index = 1: 1
    for oper_index = 1: oper_N
        lines{end+1} = ycPlot(tau_list, squeeze( OO_list(oper_index,oper_index,N_index,:)) , '$h_x\tau$', '$<O_i(t)O_i(0)>_{\beta=0}$', ...,
            ['$Dip_z+X, \frac{h_x}{J_z}=',mat2str(hx) ...
            , ',O\in \{XX,YY,ZZ\}','$'], legends, {{'line_type', line_type_list{N_index} }});
        colors{end+1} = get(lines{end}, 'color');
    end
    %legend_lines(end+1) = lines{end-3};
    %legends{end+1} = ['$N_{atom}=',mat2str(N_all(N_index)),' ,i=1','$'];
end
for N_index = 2: length(N_all)
    for oper_index = 1: oper_N
        lines{end+1} = ycPlot(tau_list, squeeze( OO_list(oper_index,oper_index,N_index,:)) + 0.02*rand(), '$h_x\tau$', '$<O_i(t)O_i(0)>_{\beta=0}$', ...,
            ['$Dip_Z+X, \frac{h_x}{J_z}=',mat2str(hx) ...
            , ',O\in \{XXX,\cdots\}','$'], legends, {{'line_type', line_type_list{N_index} }, {'color', colors{oper_index}}});
    end
    %legend_lines(end+1) = lines{end+N_index-oper_N};
    %legends{end+1} = ['$N_{atom}=',mat2str(N_all(N_index)),' ,i=', mat2str(N_index),'$'];
end
%legend(legend_lines, legends)
legend('XXX','ZXZ','YXY','XYY+YYX','XZZ+ZZX','ZYX+XYZ','XZY+YZX','ZXY+YXZ');
%ylim([0 0.2])

subplot(1,2,2)
hold on
lines = {};
for oper_index = 1: oper_N
    for N_index = 1: length(N_all) 
        lines{end+1} = ycPlot(tau_list, squeeze( OO_eig(oper_index,N_index,:))*0.4/0.568, '$h_x\tau$', 'eigs of $<O_i(t)O_j(0)>_{\beta=0}$', ...,
            ['eig values $\lambda_i$ of $O_\mu O_\nu $'], legends, {{'line_type', line_type_list{N_index} }...
            , {'color', colors{oper_index}}});
    end
end
legends = {};
for N_index = 1: length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str(N_all(N_index)),'$'];;
end
legend(legends)
%legend_lines = [lines{1} lines{1} lines{1} lines{1} lines{1} lines{1} lines{1} lines{1} ];
%legend(legend_lines, legends)
%ylim([0 0.35])
%legend('XXX','ZXZ','YXY','XYY+YYX','XZZ+ZZX','ZYX+XYZ','XZY+YZX','ZXY+YXZ');


profile viewer