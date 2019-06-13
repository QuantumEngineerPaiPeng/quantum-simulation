%20180308
% compute the ground state using exact diagonalization 
% Using H/phi to iteratively compute the ground state
% Using only matrix multiplication
TT=cputime;
N_atom=14;
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
model='dipolar';
% cplist=1;
bc='p';
glist=0.8:0.05:1.5;

H_int=-Hamiltonian(N_atom,bc,cplist,model,1);

H_trans=Sigma_z;

phi0=zeros(2^N_atom,1);
phi0(2^N_atom)=1;
% phi0=rand(2^N_atom,1);
% phi0=phi0/sqrt(phi0'*phi0);
Y2=[];
Mzlist=[];
dMzlist=[];
entropy=[];
dentropy=[];
Mylist=[];
Mxlist=[];

for g=glist
    U=-0.01*(H_int+g*H_trans)+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
    
    phi=[phi0,U*phi0];
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    
    while norm(phi(:,1)-phi(:,2))>1e-6
        phi(:,1)=phi(:,2);
        phi(:,2)=U*phi(:,2);
        phi(:,2)=phi(:,2)/norm(phi(:,2));
    end
    gs=phi(:,2);
    phi0=gs;
    Mzlist=[Mzlist,gs'*Sigma_z*gs];
    
    Ygs=Sigma_y*gs;
    Y2=[Y2,Ygs'*Ygs];
%     Gs=reshape(gs,[2^(N_atom/2),2^(N_atom/2)]);
%     s=svd(Gs);
%     entropy=[entropy,-2*log(s')*(s.^2)];
    
    U=-0.01*(H_int+(g+0.001)*H_trans)+sparse(1:2^N_atom,1:2^N_atom,ones(1,2^N_atom));
    
    phi=[phi0,U*phi0];
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    
    while norm(phi(:,1)-phi(:,2))>1e-6
        phi(:,1)=phi(:,2);
        phi(:,2)=U*phi(:,2);
        phi(:,2)=phi(:,2)/norm(phi(:,2));
    end
    gs=phi(:,2);
    phi0=gs;
    dMzlist=[dMzlist,(gs'*Sigma_z*gs-Mzlist(length(Mzlist)))/0.001];
%     Gs=reshape(gs,[2^(N_atom/2),2^(N_atom/2)]);
%     s=svd(Gs);
%     dentropy=[dentropy,(-2*log(s')*(s.^2)-entropy(length(entropy)))/0.001];
    
    g
end

cputime-TT

% figure(1)
% subplot(4,1,1)
% hold on
% plot(glist/2,Mzlist)
% 
% subplot(4,1,2)
% hold on
% plot(glist/2,entropy)
% 
% subplot(4,1,3)
% hold on
% %plot(glist(1:length(glist)-1),Mzlist(2:length(Mzlist))-Mzlist(1:length(Mzlist)-1))
% plot(glist/2,dMzlist)
% 
% subplot(4,1,4)
% hold on
% %plot(glist(1:length(glist)-1),entropy(2:length(entropy))-entropy(1:length(entropy)-1))
% plot(glist/2,dentropy)
% figure(3)
% hold on
% plot(glist,-Mzlist/N_atom)
figure(4)
hold on
plot(glist,-dMzlist/N_atom)
% figure(2)
% hold on
% plot(glist,Y2/N_atom)