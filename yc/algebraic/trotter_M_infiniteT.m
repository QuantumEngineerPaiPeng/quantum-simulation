%20190313
% compute observable ZZ(t->\infty) and XX of trotterized Ising chain with algebraically decaying
% interacions at T = \infty
% ZZ(t) = < Z(t)Z_0 >_{\beta=0} = tr( U^d(t) Z_0 U Z_0 )/ 2^N_atom, Z_0 = S^z
% H = Hz + Hx
% Hz = Jza \sum_{i<j} S^z_i S^z_j / |i-j|^alpha, 
% Jza = Jz/ (1/(N-1) * \sum_{i<j} / |i-j|^alpha )
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/16

%% runtime parameters
N_all = [8:2:12];

alpha = 3;
hx = 1;
Jz = 1;
hz = 0;

model = [1, -1/2, -1/2];
bc = 'p';
sym = ['kpz'];

dtau = 0.1;
tau_list = [dtau:dtau:3];

oper_N = 7;
t_interval = 'ED'; % only 2 element
    
profile on
profile clear

%% begin of for diff N
OO_list = deal( zeros( oper_N, length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
N_atom = N_all(N_index)
Jza = Jz; %* (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%% test operator 
Z0 = OperatorClass(N_atom,model,1/4,bc,cplist);
Z0.symmetrize(sym);
X0 = OperatorClass(N_atom,'z',1/2);
X0.symmetrize(sym);
dipx = OperatorClass(N_atom,[-1/2 ,-1/2 ,1],1/4,bc,cplist);
dipx.symmetrize(sym);
dipy = OperatorClass(N_atom,[-1/2, 1 ,-1/2],1/4,bc,cplist);
dipy.symmetrize(sym);

dipv = OperatorClass(N_atom,[0, 0, 1],1/4,bc,cplist);
dipv.symmetrize(sym);
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

opers = {Z0, X0, dipx, dipy, JP, JE, dipv};

%%  construct Hamiltonian and Floquet

U1 = H2U(Z0, dtau);
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

temp =  ZZ_longtime(Ut, t_interval, opers);
OO_list(:, N_index,tau_index) = temp*4/N_atom;

end % for tau_index = 1:length(tau_list)

end % for N_index = 1: length(N_all)

%% analytical result from PhysRevE.65.036208
analy = max( abs(cos(Jz*tau_list/2)), abs(cos(hx*tau_list)) );
analy = (analy-1)./sin(hx*tau_list).^2 +1;

%% finite size scaling
legends = cell(0);
N_all(end+1) = 16;
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end
%legends{end+1} = 'analytic';

figure(27)
clf;

subplot(1,3,1);
ycPlot(tau_list, squeeze( OO_list(1,:,:)) , '$h_x\tau$', '$<Dip_Z(n)Dip_Z(0)>_{\beta=0}$', ...,
    ['$Dip_z+X, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
subplot(1,3,2);
%}
ycPlot(tau_list,squeeze( OO_list(7,:,:)), '$h_x\tau$', '$<Dip_X(n)Dip_X(0)>_{\beta=0}$', ...,
    ['$Dip_z+X, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);

subplot(1,3,3);
%}
ycPlot(tau_list,squeeze( OO_list(4,:,:)), '$h_x\tau$', '$<Dip_Y(n)Dip_Y(0)>_{\beta=0}$', ...,
    ['$Dip_z+X, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);

%% J=0.85 vs J=1
%{
figure(92)
clf;

legends = {'$J=0.854$','$J=1$', '$J=1$ normalized'};
ycPlot(tau_list, XX_list1(end, :), '$h_x\tau$', '$<X(\infty)X(0)>_{\beta=0}$', ...,
    ['dipo, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
ycPlot(tau_list, XX_list(end-1, :), '$h_x\tau$', '$<X(\infty)X(0)>_{\beta=0}$', ...,
    ['dipo, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
ycPlot(tau_list, XX_list(end-1, :) /XX_list(end, 1) *XX_list1(end, 1) , '$h_x\tau$', '$<X(\infty)X(0)>_{\beta=0}$', ...,
    ['dipo, $N=14, \frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
%}
profile viewer