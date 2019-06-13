%20180109
% compute the ground state using exact diagonalization 
% Using H/phi to iteratively compute the ground state
TT=cputime;
N_atom=12;
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
cplist=[1];
bc='o';
glist=0.3;

H_int=sparse(2^N_atom,2^N_atom);
for pp=1:length(cplist)
    for p=1:N_atom
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else 
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
%         H_int=H_int+(-1)*cplist(pp)*...
%             (Sigma_i{1}{p}*Sigma_i{1}{q});
        q
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end

H_trans=sparse(2^N_atom,2^N_atom);
for p=1:N_atom
    H_trans=H_trans+Sigma_i{3}{p};
end

phi0=zeros(2^N_atom,1);
phi0(2^N_atom)=1;

Mzlist=[];
dMzlist=[];
entropy=[];
dentropy=[];
Mylist=[];
Mxlist=[];


ns=1;
for g=glist
    U=0.01*(H_int+g*H_trans)+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
    
    phi=[phi0,tgcr(U,phi0,tol,maxiter,phi0)];
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    
    while norm(phi(:,1)-phi(:,2))>1e-6
        phi(:,1)=phi(:,2);
        phi(:,2)=tgcr(U,phi(:,2),tol,maxiter,phi(:,2));
        phi(:,2)=phi(:,2)/norm(phi(:,2));
        ns=ns+1;
    end
    gs=phi(:,2);
    phi0=gs;
    Mzlist=[Mzlist,gs'*Sigma_z*gs];
    Gs=reshape(gs,[2^(N_atom/2),2^(N_atom/2)]);
    s=svd(Gs);
    entropy=[entropy,-2*log(s')*(s.^2)];
    
    U=0.01*(H_int+(g+0.01)*H_trans)+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
    
    phi=[phi0,tgcr(U,phi0,tol,maxiter,phi0)];
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    
    while norm(phi(:,1)-phi(:,2))>1e-6
        phi(:,1)=phi(:,2);
        phi(:,2)=tgcr(U,phi(:,2),tol,maxiter,phi(:,2));
        phi(:,2)=phi(:,2)/norm(phi(:,2));
    end
    gs=phi(:,2);
    phi0=gs;
    dMzlist=[dMzlist,(gs'*Sigma_z*gs-Mzlist(length(Mzlist)))/0.01];
    Gs=reshape(gs,[2^(N_atom/2),2^(N_atom/2)]);
    s=svd(Gs);
    dentropy=[dentropy,(-2*log(s')*(s.^2)-entropy(length(entropy)))/0.01];
    
    g
end

cputime-TT

figure(1)
subplot(4,1,1)
hold on
plot(glist/2,Mzlist)

subplot(4,1,2)
hold on
plot(glist/2,entropy)

subplot(4,1,3)
hold on
%plot(glist(1:length(glist)-1),Mzlist(2:length(Mzlist))-Mzlist(1:length(Mzlist)-1))
plot(glist/2,dMzlist)

subplot(4,1,4)
hold on
%plot(glist(1:length(glist)-1),entropy(2:length(entropy))-entropy(1:length(entropy)-1))
plot(glist/2,dentropy)