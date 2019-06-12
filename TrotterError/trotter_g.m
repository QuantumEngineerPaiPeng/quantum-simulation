%20180711
% compute infinite time averaged operator for different floquet unitary 
% Written using OperatorClass
% Created by Pai Peng
N_atom=12;
tic

bc='p';
% cplist=[1,2,3].^-3;
cplist=1;
% cplist=(1:N_atom-1).^-3;
% model='Ising';
model='dipolar';
% sym='kpz';
sym='kp';
% sym=[];
% tlist=(0:0.6281:10.5)/4;
% tlist=logspace(-1,5,50);
% tlist=linspace(0,5,100);
% glist=linspace(2,10,5); % strength of the transverse field
% glist=[3];
% glist=(0:0.2:2);
glist=logspace(-1,2,30);
T=2*pi/0.2;

rho_0=OperatorClass(N_atom,'z',1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,1,'x');
V0=(copy(rho_0));

rho_0.symmetrize(sym);
V0.symmetrize(sym);

C=zeros(length(glist),1);

H_int=OperatorClass(N_atom,model,-1,bc,cplist,1);
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1);
Hx.symmetrize(sym);

H_int=H_int+2*Hz;

col=1;
parfor p=1:length(glist)
    g=glist(p);
    H1=H_int+g*Hz;
    H2=H_int-g*Hz;

    U=H2U(H1,T/2)*H2U(H2,T/2);
    H=U2H(U,T);
%     C(:,p)=TPC(rho_0,H,V0,tlist);
%     C(:,p)=OTOC(rho_0,H,V0,tlist);
%     C(:,p)=BiEETlist(rho_0,H,tlist,5);
    rho=InfAveOpe(rho_0,H);
    C(p)=overlap(rho,rho);
end

C=C/(N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure(6)
plot(glist,real(C))
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