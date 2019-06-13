% 20180510
% Using the saved xx prethermal Hamiltonian to calculate r as a function of
% epsilon
function prexx_r_order_func(N_atom,order,epslist)
% tic
% N_atom=8;
% J=0.25;
bc='o';
fname=...
    '~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/prethermalHxx';

load([fname,'/hpre_xx',num2str(N_atom)],'hpre')
if N_atom==11
load([fname,'/2018_5_10_9_38_7_.mat'])
else if N_atom==8
        load([fname,'/2018_5_9_20_55_15_.mat'])
    end
end
% order=length(hpre)-1;
% order=6;
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
% epslist=logspace(-1,1,10);
% epslist=0.25:0.25:3;
r=zeros(order+1,length(epslist));

pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
E0=zeros(length(epslist),2^N_atom);
for p=1:order+1
    rtemp=zeros(1,length(epslist));
    parfor count=1:length(epslist)
        eps=epslist(count); % eps is the prefactor of interaction in spin picture
        H=H0+eps*H1;
        if p==1
            E0temp=sort(eig(full(H)));
            E0(count,:)=E0temp;
        end
        Hpre=zeros(2^N_atom,2^N_atom);
        for pp=1:p
            Hpre=Hpre+eps^(pp-1)*hpre{pp};
        end
        rtemp(count)=mean(abs(E0(count,:)-sort(real(eig(full(Hpre)))).'))/N_atom;
    end
    r(p,:)=rtemp;
    save([filetosave,'r_',num2str(N_atom),num2str(order)],'r','N_atom','p','epslist','cplist');
end
% toc


figure(1)
semilogy(0:order,r)
xlabel('order')
ylabel('r')
ppStyle(20,2)