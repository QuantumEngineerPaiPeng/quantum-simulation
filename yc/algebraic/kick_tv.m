%20190531
% repeat PRE 60, 3949, with the view of point of trotterized Hessenberg
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/31

%% runtime parameters
N_all = [10];

alpha = 100;
hx = 1; % V/t
Jz = 1;

bc = 'p';
sym = [];

dtau = 1;
tau_list = [0:dtau:8];
n_end = 100;
dn = 5;
n_list = [1: dn: n_end];

%F_N = [10^2];

load_save = false;
diago = true;
    
profile on
profile clear

%% begin of for diff N
[ZZ_list, XX_list, JJ_list] = deal( zeros(length(n_list), length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
N_atom = N_all(N_index)
Jza = Jz; %* (N_atom-1) / ( sum( [N_atom-1:-1:1] ./ [1:N_atom-1].^alpha  ) );
cplist = Jza * [1:N_atom-1].^(-alpha) ;

%% test operator 
Z0 = OperatorClass(N_atom,[1, 1, 0],1/4,bc,cplist);
X0 = OperatorClass(N_atom,[0, 0, 1],1/4,bc,cplist) ;% + OperatorClass(N_atom,'z',1/2); this term is irrelevant
Z0.symmetrize(sym);
X0.symmetrize(sym);

%% test operator: current J = XY-YX (when dipo is Z), here = YZ-ZY
J0 = OperatorClass(N_atom);
J0.matrix = { sparse(2^N_atom, 2^N_atom) };
for atom_index = 1: N_atom
    J0.matrix{1} = J0.matrix{1} + Kron2body(N_atom,1, 2,atom_index, mod(atom_index,N_atom)+1)/2 ...
        - Kron2body(N_atom,2, 1,atom_index, mod(atom_index,N_atom)+1)/2;
end
J0.symmetrize(sym);

%%  construct Hamiltonian and Floquet
Hz= Z0; %+ hz* OperatorClass(N_atom,'x',1/2);
%Hz.symmetrize(sym);

Hx=hx*	X0;
%Hx.symmetrize(sym);

U1 = H2U(Hz, dtau);
U2 = H2U(Hx, dtau);

%% begin of for diff tau
for tau_index = 1:length(tau_list)

tau = tau_list(tau_index)
if tau < 0.001
    Ut = H2U( Hz+ Hx, 1);
elseif tau_index == 1
    Ut = U1* U2;
else
    Ut = U1 * Ut * U2;
end

temp =  ZZ_longtime(Ut, [1, n_end], {Z0, X0, J0}, {{'all_t',true}, {'dn', dn}});
ZZ_list(:, N_index,tau_index) = temp(1, :)*4/N_atom;
XX_list(:, N_index,tau_index) = temp(2, :)*4/N_atom;
JJ_list(:, N_index,tau_index) = temp(3, :)*4/N_atom /2;

end % for tau_index = 1:length(tau_list)

end % for N_index = 1: length(N_all)

%% C(J) vs n
legends = cell(0);
for N_index = 1:length(tau_list)
    legends{end+1} = ['$t\tau=',mat2str( tau_list(N_index) ),'$'];
end
figure(25)
clf;

subplot(1,2,1);
ycPlot(n_list, squeeze( ZZ_list(:,end,:)), '$n_{period}$', '$<H_1(n)H_1(0)>_{\beta=0}$', ...,
    ['kick-tV, $\frac{V}{t}=',mat2str(hx) ...
    , '$'], legends);
subplot(1,2,2);
ycPlot(n_list, squeeze( JJ_list(:,end,:)), '$n_{period}$', '$<J(n)J(0)>_{\beta=0}$', ...,
    ['kick-tV, $\frac{V}{t}=',mat2str(hx) ...
    , '$'], legends);

%% finite size scaling
%{
legends = cell(0);
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end
%legends{end+1} = 'analytic';

figure(27)
clf;

subplot(1,2,1);
ycPlot(tau_list, ZZ_list, '$h_x\tau$', '$<H_{dipz}(\infty)H_{dipz}(0)>_{\beta=0}$', ...,
    ['dipo, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
subplot(1,2,2);

ycPlot(tau_list, JJ_list, '$h_x\tau$', '$<X(\infty)X(0)>_{\beta=0}$', ...,
    ['dipo, $\frac{h_x}{J_z}=',mat2str(hx),' ,\frac{h_z}{J_z}=',mat2str(hz) ...
    , '$'], legends);
%}
profile viewer