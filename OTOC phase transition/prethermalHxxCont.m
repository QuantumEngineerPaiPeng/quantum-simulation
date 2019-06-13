% 20180509
% compute the prethermal Hamiltonian from saved s and hpre matrices
% function prethermalHxx_func(N_atom,order,bc,cplist,g)
tic
N_atom=12;
g=0;
J=0.25;
bc='o';
order=N_atom;
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
clear('Sigma_i')
cplist=(1:N_atom-1).^(-3);
% cplist=1;

H1=-J*Hamiltonian(N_atom,bc,cplist,'dipolar',1)+g*Sigma_z;

H0=-J*Hamiltonian(N_atom,'o',1,'Ising',1);

H1=H1-H0;
% Hpre=zeros(2^N_atom,2^N_atom,order);
% Hpre=H0;
% hpre={H0};
% s={};
% s_sum=zeros(2^N_atom,2^N_atom);

% vanM=zeros(2*N_atom,2*N_atom);
% temp=[-N_atom:-1,1:N_atom];
% for p=1:2*N_atom
%     vanM(p,:)=temp.^p;
% end
% ivanM=inv(vanM);
filetoload=...
    '~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/prethermalHxx/';
hprepath=[filetoload,'2018_5_9_9_44_21_hpre_xx12.mat'];
spath=[filetoload,'2018_5_9_9_44_21_s_xx12.mat'];
load(hprepath)
load(spath)
load('~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ivanM.mat',['ivanM',num2str(N_atom)])
ivanM=eval(['ivanM',num2str(N_atom)]);
% r=zeros(1,order+1);
% rH=r;
% E0=sort(eig(H0+H1));
% r(1)=mean(abs(E0-sort(real(eig(H0)))))/N_atom;
% rH(1)=sum(sum(abs(H1.^2)));

for p=length(s)+1:order
    h=expandSH(s,H0,p)+expandSH(s,H1,p-1); % Hpre(n)=[sp,H0]+h
    zh={}; % list contains [H0,h],[H0,[H0,h]],...
    temp=h;
    for pp=1:2*N_atom
        temp=(H0*temp-temp*H0)/J*0.5;
        zh{end+1}=temp;
    end
    mqc={};
    for pp=1:2*N_atom
        temp=sparse(2^N_atom,2^N_atom);
        for ppp=1:2*N_atom
            temp=temp+ivanM(pp,ppp)*zh{ppp};
        end
        mqc{end+1}=temp;
    end
    
    temp=sparse(2^N_atom,2^N_atom);
    for pp=1:length(mqc)
        temp=temp+mqc{pp};
    end
%     Hpre=Hpre+h-temp;
    hpre{end+1}=h-temp;
    if g==0
    save(hprepath,'hpre','p')
    end
%     r(p+1)=mean(abs(E0-sort(real(eig(full(Hpre))))))/N_atom;
    
    count=1;
    temp=sparse(2^N_atom,2^N_atom);
    for qc=[-N_atom:-1,1:N_atom]
        temp=temp+mqc{count}/qc;
        count=count+1;
    end
%     s_sum=s_sum+temp;
    s{end+1}=temp;
    save(spath,'s','p')
%     rH(p+1)=sum(sum(abs((expm(s_sum)*(H0+H1)*expm(-s_sum)-Hpre).^2)));
%     save(filetosave,'r','rH','N_atom','g','J','bc','cplist')
end
toc
