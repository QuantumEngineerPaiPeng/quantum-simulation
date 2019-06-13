%20190522
% compute auto-correlation in PhysRevX.4.041048
% Written using OperatorClass
% Created by Chao Yin
% last change: 5/22

%% runtime parameters
N_all = [8:2:14];

J = 1;
dJ = 0.2;
Jp = 0.8;

bc = 'p';
sym = ['kpz'];

dtau = 0.05;
tau_list = [dtau:dtau:2];

F_N = [10^2];

load_save = false;
diago = true;
    
profile on
profile clear

%% begin of for diff N
[ZZ_list, XX_list] = deal( zeros(length(N_all), length(tau_list)) );
for N_index = 1: length(N_all)
    
N_atom = N_all(N_index)

%%  construct Hamiltonian and Floquet
dipz = OperatorClass(N_atom, [-1/2, -1/2 ,1] ,1/4,bc,[1]);
dipz.symmetrize(sym);

ising_next = OperatorClass(N_atom, [0, 0 ,1] ,1/4,bc,[0 2]);
ising_next.symmetrize(sym);

H1 = (J+dJ)* dipz + Jp* ising_next;
H2 = H1 - 2*dJ* dipz;

if load_save && exist(filename,'file')
    load(filename);
    loaded = true
else
    Ut_all = cell(0);
    Ut=OperatorClass(N_atom);
    Ut.matrix = { speye(2^N_atom) };
    Ut.symmetrize(sym);
    U1 = H2U(H1, dtau);
    U2 = H2U(H2, dtau);
    
    loaded = false
end % if exist(filename,'file')

%% test operator 
Z0 = H1;% OperatorClass(N_atom,'z',1/2);
X0 = dipz;
%Z0.symmetrize(sym);
%X0.symmetrize(sym);

%% begin of for diff tau
for tau_index = 1:length(tau_list)

tau = tau_list(tau_index);
if loaded
    Ut = Ut_all{tau_index};
else
    Ut = U1 * Ut * U2;
    %Ut_all{end+1} = copy( Ut);
end

if diago
    F_N = 'ED';
end
temp =  ZZ_longtime(Ut, F_N, {Z0,X0});
ZZ_list(N_index,tau_index) = temp(1, :)*4/N_atom;
XX_list(N_index,tau_index) = temp(2, :)*4/N_atom;

end % for tau_index = 1:length(tau_list)

%% after for, save
if load_save && ~loaded
    save(filename, 'Ut_all');
end

end % for N_index = 1: length(N_all)


%% finite size scaling
legends = cell(0);
for N_index = 1:length(N_all)
    legends{end+1} = ['$N_{atom}=',mat2str( N_all(N_index) ),'$'];
end

figure(26)
clf;

subplot(1,2,1);
ycPlot(tau_list, ZZ_list, '$J_z\tau$', '$<H_1(\infty)H_1(0)>_{\beta=0}$', ...,
    ['Rigol, $J,\delta J,\tilde{J}=',mat2str(J),',',mat2str(dJ),',',mat2str(Jp) ...
     ' ,N_{period}=\infty$'], legends);
subplot(1,2,2);
%}
ycPlot(tau_list, XX_list, '$h_x\tau$', '$<Dipz(\infty)Dipz(0)>_{\beta=0}$', ...,
    ['Rigol, $J,\delta J,\tilde{J}=',mat2str(J),',',mat2str(dJ),',',mat2str(Jp) ...
     ' ,N_{period}=\infty$'], legends);

profile viewer