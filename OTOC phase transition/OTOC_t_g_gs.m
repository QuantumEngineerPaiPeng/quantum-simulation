%20180114
% compute the ground state using exact diagonalization
% Using H/phi to iteratively compute the ground state
% then compute OTOC at the ground state
TT=cputime;
N_atom=11;
N_rand=1;
tol=1e-8;
maxiter=100;
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end

glist=logspace(-1,1,10);
tlist=logspace(2,5,20);
C=zeros(length(tlist),length(glist));
rho_0=Sigma_y;
V0=Sigma_z;

H_int=sparse(2^N_atom,2^N_atom);

for p=1:N_atom-1
    H_int=H_int+(-1)*...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
%         H_int=H_int+(-1)*...
%             (Sigma_i{1}{p}*Sigma_i{1}{p+1});
    %     H_int=H_int+...
    %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
    
end

H_trans=sparse(2^N_atom,2^N_atom);
for p=1:N_atom
    H_trans=H_trans+Sigma_i{3}{p};
end

phi0=zeros(2^N_atom,1);
phi0(2^N_atom)=1;

Mzlist=[];

col=1;
for g=glist
    U=0.001*(H_int+g*H_trans)+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
    
    phi=[phi0,tgcr(U,phi0,tol,maxiter)];
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    
    while norm(phi(:,1)-phi(:,2))>1e-6
        phi(:,1)=phi(:,2);
        phi(:,2)=tgcr(U,phi(:,2),tol,maxiter);
        phi(:,2)=phi(:,2)/norm(phi(:,2));
    end
    gs=phi(:,2);
    phi0=gs;
    Mzlist=[Mzlist,gs'*Sigma_z*gs];
    H=full(H_int+g*H_trans);
    [V,D]=eig(H);
    row=1;
    for t=tlist
        rhot=V*diag(exp(diag(1i*D*t)))*V'*rho_0*V*diag(exp(diag(-1i*D*t)))*V';
        C(row,col)=-gs'*((rhot*V0-V0*rhot)^2)*gs;
        row=row+1;
    end
    col=col+1;
    g
end

cputime-TT


figure(3)
hold on
plot(glist,C/(2*N_atom*2^N_atom))
