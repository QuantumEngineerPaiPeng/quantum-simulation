%20180308
% compute the ground state using exact diagonalization
% Using H/phi to iteratively compute the ground state
% Using only matrix multiplication
TT=cputime;
N_atom=8;
N_rand=1;
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
cplist=[1,1/8,1/27];
% cplist=1;
model='dipolar';
bc='p';
glist=0.02:0.2:3;
tlist=0.:0.5:3;
clear('Sigma_i')

H_int=-Hamiltonian(N_atom,bc,cplist,model,1);
V1=Sigma_z;
V2=Sigma_x;
phi0=zeros(2^N_atom,2);
phi0(2^N_atom,1)=1;
phi0(2^N_atom-1,2)=1;
Y2=[];
E=zeros(length(glist),2);
C=zeros(length(tlist),length(glist));
for col=1:length(glist)
    for sec=1:2
        H=H_int+glist(col)*Sigma_z;
        U=-0.01*H+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
        phi=[phi0(:,sec),U*phi0(:,sec)];
        phi(:,2)=phi(:,2)/norm(phi(:,2));
        while norm(phi(:,1)-phi(:,2))>1e-8
            phi(:,1)=phi(:,2);
            phi(:,2)=U*phi(:,2);
            phi(:,2)=phi(:,2)/norm(phi(:,2));
        end
        gs=phi(:,2);
        phi0(:,sec)=gs;
        E(col,sec)=gs'*H*gs;
    end
    rho=1/2*phi(:,1)*phi(:,1)'+1/2*phi(:,2)*phi(:,2)';
    [V,D]=eig(full(H));
    for row=1:length(tlist)
        t=tlist(row);
        V1t=V*diag(exp(diag(1i*D*t)))*V'*V1*V*diag(exp(diag(-1i*D*t)))*V';
        C(row,col)=trace(rho*(V1t*V2-V2*V1t)*(V1t*V2-V2*V1t)');
    end
end
C=C/(2^N_atom);
figure(1)
hold on
plot(glist,C)
cputime-TT