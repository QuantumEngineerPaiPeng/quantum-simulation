%20180620
% Magnetizatino for different t
% Using symmetries to speed up
% Written using operator class
% Evolution calculated using the power of Floquet Unitary
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

H_int=OperatorClass(L, model, -1, bc, cplist, 1);
H_int.symmetrize(sym);
H_int.diagonalize;

% toc
% filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
% if exist(filename,'file')
%     load(filename)
% else
%     pt=AllkpzProject(N_atom);
%     save(filename,'pt')
% end
% toc

rho0=OperatorClass(L, 'x', 1);
rho0.symmetrize(sym);

Hp=OperatorClass(L, 'z', 1);
Hp.symmetrize(sym);
Hp.diagonalize;

M=zeros(length(epslist),length(Nlist));

for p=1:length(epslist)
    eps=epslist(p);
    Up=H2U(Hp,eps/2);
    UH=H2U(H_int,T);
    Uf=Up*UH;
    temp=zeros(length(Nlist),1);
    for q=1:length(Nlist)
        N=Nlist(q);
        UN=Uf^N;
        temp(q)=overlap(UN*rho0*UN',rho0);
    end
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
plot(Nlist,Mt,'LineWidth',2)
xlabel('# of Floquet periods')
ylabel('|M_z|')
