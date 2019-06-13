%20180620
% Magnetizatino for different t
% Using symmetries to speed up
% Written using operator class
% Evolution calculated using Floquet Hamiltonian
% created by Pai Peng
% function M_t_sym_func2(L,Nlist,epslist,T,bc,cplist,model,sym)
tic
L=12;
Nlist=1:1:100;
epslist=0.1;
model='dipolar';
T=0.5;
bc='p';
cplist=[1,1/8,1/27];
sym='kp';

% H_int=-Hamiltonian(N_atom,bc,cplist,model,1);

H_int=OperatorClass(L, model, -1, bc, cplist, 1);
H_int.symmetrize(sym);
H_int.diagonalize;

rho_0=OperatorClass(L, 'x', 1);
rho_0.symmetrize(sym);

Hp=OperatorClass(L, 'z', 1);
Hp.symmetrize(sym);
Hp.diagonalize;

M=zeros(length(epslist),length(Nlist));

for p=1:length(epslist)
    eps=epslist(p);
    Up=H2U(Hp,eps/2);
    UH=H2U(H_int,T);
    Uf=UH*Up;
    Hf=U2H(Uf,T);
    temp=TPC(rho_0,Hf,rho_0,T*Nlist);
    M(p,:)=temp;
end



Mt=M/(L*2^L);
toc

% pfolder='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
% foldername=mfilename;
% mkdir(pfolder,foldername)
% filetosave=[pfolder,foldername,'/'];
% for p = fix(clock)
%     filetosave=[filetosave,num2str(p),'_'];
% end
% filetosave=[filetosave,'.mat'];
% save(filetosave,'N_atom','model','epslist','Nlist','T','bc','cplist','Mt')
% varname={'file','N_atom','model','epslist','Nlist','T','bc','cplist'};
% varvalue={filetosave,N_atom,model,epslist,Nlist,T,bc,cplist};
% writelog(pfolder,varname,varvalue)
figure(2)
hold on
plot(Nlist,Mt)
xlabel('# of Floquet periods')
ylabel('|M_z|')
