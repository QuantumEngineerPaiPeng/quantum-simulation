% 20190701
% Calculate the FID dipolar Z state
% created by Pai Peng
N_atom=10;
tic

bc='p';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
% model='Ising';
model=[-0.5,-0.5,1];
% sym='kpz';
sym='kp';
% sym=[];
% tlist=(0:0.6281:10.5)/4;
% tlist=logspace(-1,5,50);
tlist=0:0.05:3;
q_max=8;

rho_0=OperatorClass(N_atom,model,-1,bc,cplist);
rho_0.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1);
Hx.symmetrize(sym);
U45=H2U(Hx,pi/4/2);
rho_0=U45*rho_0*U45';

V0=OperatorClass(N_atom,'x',1);
V0.symmetrize(sym);

% C=zeros(length(tlist),1);

H_int=OperatorClass(N_atom,model,-1,bc,cplist);
H_int.symmetrize(sym);

% col=1;

C=TPC(rho_0,H_int,V0,tlist);

% C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure
plot(tlist,real(C)/(N_atom*2^N_atom))
hold on
% xlabel('g')
ppStyle(30,2,10)

% pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
% foldername=mfilename;
% mkdir(pfolder,foldername)
% filetosave=[pfolder,foldername,'/'];
% for p = fix(clock)
%     filetosave=[filetosave,num2str(p),'_'];
% end
% filetosave=[filetosave,'.mat'];
% save(filetosave,'N_atom','model','glist','tlist','rho_0code','V0code','bc','cplist','C')
% varname={'file','N_atom','model','glist','tlist','rho_0code','V0code','bc','cplist'};
% varvalue={filetosave,N_atom,model,glist,tlist,rho_0code,V0code,bc,cplist};
% writelog(pfolder,varname,varvalue)
toc