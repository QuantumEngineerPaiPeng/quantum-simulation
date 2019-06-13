% 20180505
% compute the prethermal Hamiltonian
function prethermalH_func(N_atom,order,bc,cplist)
% N_atom=9;
J=0.25;
% bc='o';
% order=N_atom;
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
% cplist=(1:N_atom-1).^(-3);
% cplist=1;

H1=-J*Hamiltonian(N_atom,bc,cplist,'dipolar',1);

H0=0.5*Sigma_z;

% Hpre=zeros(2^N_atom,2^N_atom,order);
Hpre=H0;
hpre={H0};
s={};

% vanM=zeros(2*N_atom,2*N_atom);
% temp=[-N_atom:-1,1:N_atom];
% for p=1:2*N_atom
%     vanM(p,:)=temp.^p;
% end
% ivanM=inv(vanM);
load('~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ivanM.mat',['ivanM',num2str(N_atom)])
ivanM=eval(['ivanM',num2str(N_atom)]);
r=zeros(1,order+1);
rH=r;
E0=sort(eig(H0+H1));
r(1)=sum(abs(E0-sort(eig(H0))));
rH(1)=sum(sum(abs(H1.^2)));

pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];

for p=1:order
    h=expandSH(s,H0,p)+expandSH(s,H1,p-1); % Hpre(n)=[sp,H0]+h
    zh={}; % list contains [H0,h],[H0,[H0,h]],...
    temp=h;
    for pp=1:2*N_atom
        temp=(H0*temp-temp*H0);
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
    Hpre=Hpre+h-temp;
    hpre{end+1}=h-temp;
    save([filetosave,'hpre',num2str(N_atom)],'hpre','p')
    r(p)=sum(abs(E0-sort(eig(full(Hpre)))));
    
    count=1;
    temp=sparse(2^N_atom,2^N_atom);
    for qc=[-N_atom:-1,1:N_atom]
        temp=temp+mqc{count}/qc;
        count=count+1;
    end
    s{end+1}=temp;
    save([filetosave,'s',num2str(N_atom)],'s','p')
%     rH(p)=sum(sum(abs((expm(sum(s,3))*(H0+H1)*expm(-sum(s,3))-sum(Hpre,3)).^2)));
end
% toc
% figure(1)
% hold on
% semilogy(0:order,r)
% figure(2)
% semilogy(0:order,rH)