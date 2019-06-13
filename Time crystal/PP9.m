%20180624
% Magnetization for different t. Simulate PP9 sequence.
% Using symmetries to speed up
% Written using operator class
% Evolution calculated using Floquet Hamiltonian
% created by Pai Peng
% function M_t_sym_func2(L,Nlist,epslist,T,bc,cplist,model,sym)
tic

L=10;
Nlist=0:2:256;
epslist=;
model='dipolar';
cycle=36;
J=8.18e-3;
p2=2;
bc='p';
cplist=[1,1/8,1/27];
% cplist=1;
sym='kp';

tau=cycle/9-p2;

H_int=OperatorClass(L, model, -J, bc, cplist, 3);
H_int.symmetrize(sym);
H_int.diagonalize;

rho_0=OperatorClass(L, 'z', 1);
rho_0.symmetrize(sym);

phase=[0,60,0,120,180,120,0,60,0];
delay=tau*[0.5,1,1,1,1,1,1,1,1,0.5];

M=zeros(length(epslist),length(Nlist));

for p=1:length(epslist)
    eps=epslist(p);
    angle=[pi*ones(1,8),pi+eps];
    Uf=UFloquet(phase,angle,delay,H_int,p2);
    Hf=U2H(Uf,1);
    
    M(p,:)=TPC(rho_0,Hf,rho_0,Nlist);
end
Mt=M/(L*2^L);
toc

pfolder='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
for p = fix(clock)
    filetosave=[filetosave,num2str(p),'_'];
end
filetosave=[filetosave,'.mat'];
save(filetosave,'L','model','epslist','Nlist','cycle','J','p2','bc','cplist','sym','Mt')
varname={'filetosave','L','model','epslist','Nlist','cycle','J','p2','bc','cplist','sym'};
varvalue=cell(size(varname));
for p=1:length(varname)
    varvalue{p}=eval(varname{p});
end
writelog(pfolder,varname,varvalue)
% figure(1)
% hold on
% plot(Nlist,Mt)
% xlabel('# of Floquet periods')
% ylabel('|M_z|')
