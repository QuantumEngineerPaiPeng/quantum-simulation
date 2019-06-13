%20180622
% Magnetization for different t. Simulate master sequence.
% Using symmetries to speed up
% Written using operator class
% Evolution calculated using Floquet Hamiltonian
% created by Pai Peng
% function M_t_sym_func2(L,Nlist,epslist,T,bc,cplist,model,sym)
tic

L=8;
Nlist=0:1:100;
eps=0.;
model='dipolar';
cycle=100;
J=8.18e-3;
p1=1;
bc='p';
% cplist=[1,1/8,1/27];
cplist=1;
% sym='kp';
sym=[];

tau=cycle/24;
u=-0.2;
v=0.2;
w=0;

t1=tau*(1-v+w)-p1/2;
t2=tau*(1-u+v)-p1;
t3=tau*(1+u-w)-p1/2;

H_int=OperatorClass(L, model, -J, bc, cplist, 3);
H_int.symmetrize(sym);
H_int.diagonalize;

rho_0=OperatorClass(L, 'y', 1);
rho_0.symmetrize(sym);

phase=[0,90,90,0,0,90,90,0,180,270,270,180,180,270,270,180];
angle=pi/2;
delay=[t1,t2,2*t3,t2,2*t1,t2,2*t3,t2,2*t1,t2,2*t3,t2,2*t1,t2,2*t3,t2,t1];
Z=OperatorClass(L,'z',1);
Z.symmetrize(sym);
Uf=H2U(Z,pi*(1+eps)/2)*UFloquet(phase,angle,delay,H_int,p1);
Hf=U2H(Uf,1);

M=TPC(rho_0,Hf,rho_0,Nlist);

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
figure(1)
hold on
plot(Nlist,abs(Mt))
xlabel('# of Floquet periods')
ylabel('|M_z|')
