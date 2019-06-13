%20180108
%compute 1-body term of O(t)
tic
N_atom=13;
N_rand=1;
cplist=[1,1/8,1/27,1/4^3,1/5^3];
% cplist=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

tlist=logspace(-1,3,40);

glist=[0.25,0.5,1,2];
dg=1e-1;
gdg=reshape([glist;glist+dg],1,2*length(glist));
% rho_0=Sigma_i{3}{round((N_atom+1)/2)};
rho_0=Sigma_x;
% rho_0=Hamiltonian(N_atom,'o',1,'Ising',1);

%RandomM=rand(N_rand,3,N_atom);
bc='o';
% f1=zeros(length(gdg),length(tlist));
f1=zeros(length(glist),length(tlist));
H_int=-Hamiltonian(N_atom,bc,cplist,'dipolar',1);

row=1;
% for g=gdg
pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
for p = fix(clock)
    filetosave=[filetosave,num2str(p),'_'];
end

for g=glist
    H_trans=g*Sigma_z;
    H=H_int+H_trans;
    [V,D]=eig(H);
    col=1;
    for t=tlist
        rho=V*diag(exp(diag(1i*D*t)))*V'*rho_0*V*diag(exp(diag(-1i*D*t)))*V';
        for p=1:N_atom
            A=rho(1:2^(N_atom-p),1:2^(N_atom-p));
            B=rho(2^(N_atom-p)+1:end,1:2^(N_atom-p));
            D1=rho((1+2^(N_atom-p)):end,(1+2^(N_atom-p)):end);
            f1(row,col)=f1(row,col)+abs(2*trace(B))^2+abs(trace(A-D1))^2;
            rho=A+D1;
        end
        
        col=col+1;
        
    end
    row=row+1;
    g
end

f1=f1/(N_atom*4^N_atom);
save([filetosave,'.mat'],'N_atom','glist','tlist','bc','cplist','f1','rho_0')


% figure(1)
% semilogx(tlist,f1)
% f1_ave=real(sum(f1,2))'/length(tlist);
% df1=(f1(2:2:end,:)-f1(1:2:end,:))/dg;
% errf1=std(f1,0,2)/sqrt(length(tlist));
% errdf1=std(df1,0,2)/sqrt(length(tlist));
% figure(1)
% errorbar(gdg,f1_ave,errf1/2,'LineWidth',2)
% xlabel('g')
% ylabel('f1')
% hold on
% figure(2)
% errorbar(glist+dg/2,mean(df1,2),errdf1/2,'LineWidth',2)
% xlabel('g')
% ylabel('df1/dg')
% hold on
toc