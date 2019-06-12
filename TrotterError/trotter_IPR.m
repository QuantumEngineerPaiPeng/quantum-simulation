%20190302
% compute IPR of trotterized traverse Ising chain
% H = HZ + Hx
% HZ = \sum S^z_j S^z_{j+1} + h \sum S^z_j  (J=1)
% Hx = g \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin

%% runtime parameters
N_atom = 14;

model = [0, 0, 1];
bc = 'p';
cplist = [1];
sym = ['kp'];

g = 2;
h = 2;

IPR_ED = true;
IPR_N = 10^4;

dtau = 0.1;
tau_list = [dtau:dtau:2];

%% test state of IPR
% psi0 = zeros(2^N_atom ,1);
% psi0(1) = 1;
% psi0_t = psi0';
psi0=OperatorClass(N_atom);
psi0matrix=sparse(2^N_atom,2^N_atom);
psi0matrix(1,1)=1;
psi0.matrix={psi0matrix};
psi0.symmetrize(sym);

%% construct Hamiltonian before trotter
H_int=OperatorClass(N_atom,model,1/4,bc,cplist);  % H_int = \sum \sigma^z_j \sigma^z_{j+1}
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1/2);  % Hz = \sum \sigma^z_j
Hz.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1/2);
Hx.symmetrize(sym);

HZ = H_int + h*Hz;

H = HZ + g*Hx;

% Ut=OperatorClass(N_atom);
% Ut.matrix = { eye(2^N_atom) };

%% U1, U2 for dtau
U1 = H2U( HZ, dtau );
U2 = H2U( g * Hx, dtau );

Ut=U1*U2;

%% begin of for
IPR_list = tau_list;

for tau_index = 1:length(tau_list)
    
    tau = tau_list(tau_index)
    
    
    
    %% IPR of psi0
    tic;
    if IPR_ED == true
        Ut.diagonalize();
        temp=0;
        for p=1:length(Ut.matrix)
            for pp=1:size(Ut.matrix{p},2)
                psi_nv=Ut.eigsys{p}.V(:,pp);
                temp=temp+(psi_nv'*psi0.matrix{p}*psi_nv)^2;
            end
        end
        IPR_list(tau_index)=temp;
        if tau_index < length(tau_list)
            Ut = U1 * Ut * U2;
        end
        %     IPR_list(tau_index) = sum( abs( psi0_t * Ut.eigsys{1}.V).^4 );
    else
        psi = Ut.matrix{1}^IPR_N * psi0;
        temp = 0;
        for t = 1:IPR_N
            psi = Ut.matrix{1} * psi;
            temp = temp + abs( psi0_t * psi )^2;
        end % for t = 1:IPR_N
        
        IPR_list(tau_index) = temp / IPR_N;
    end
    
    
end % for len_index = 1:length(N_atom_list)

%% after for
log_IPR = -log2(IPR_list) ./ ( N_atom - 1 )

%% plot log_IPR vs N_atom
% ycPlot(2, tau_list, log_IPR, '$\tau$', '$\lambda_{IPR} / \lambda_{D}$', ...,
%     ['$g=',mat2str(g),' ,h=',mat2str(h),' ,N_{atom}=',mat2str(N_atom),'$']);
figure
plot(tau_list,log_IPR)