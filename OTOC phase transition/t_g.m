%20180711
% compute two-point correlator or OTOC with trasverse field for different field strength
% and time
% Written using OperatorClass
% Created by Pai Peng
N_atom=10;
tic

bc='p';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
% model='Ising';
model=[1,-1,0]*0+[-0.5,1,-0.5];
% sym='kpz';
sym='kp';
% sym=[];
% tlist=(0:0.6281:10.5)/4;
% tlist=logspace(-1,5,50);
tlist=0:1:320;
% glist=linspace(2,10,5); % strength of the transverse field
% glist=[3];
glist=[2];
% glist=logspace(-1,2,10);
theta=0/180*pi;

rho_0=OperatorClass(N_atom,'z',1)*cos(theta)+OperatorClass(N_atom,'x',1)*sin(theta);
% rho_0.matrix={zeros(2^N_atom,2^N_atom)};
% dir1=[90,0];
% dir2=[90,90];
% rho_0=StateN({dir1,dir2,dir1,dir2,dir1,dir2,dir1,dir2,dir1,dir2});
% rho_0.matrix{1}(1,1)=1;
% Hy=OperatorClass(N_atom,'y',1);
% rho_0=H2U(Hy,pi/8)*OperatorClass(N_atom,model,-1,bc,cplist,3)*(H2U(Hy,pi/8))';
% rho_0=OperatorClass(N_atom,model,-1,bc,cplist,2);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,1,'x');
V0=OperatorClass(N_atom,'z',1)*cos(theta)+OperatorClass(N_atom,'x',1)*sin(theta);

rho_0.symmetrize(sym);
V0.symmetrize(sym);

C=zeros(length(tlist),length(glist));

H_int=OperatorClass(N_atom,model,-1,bc,cplist,2);
H_int.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1);
Hx.symmetrize(sym);
Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);
col=1;
parfor p=1:length(glist)
    g=glist(p);
    H=H_int+g*Hz;
%     H=U2H(H2U(0.25*H_int,3/180*pi)*H2U(0.5*Hz,3/180*pi),1);
    C(:,p)=TPC(rho_0,H,V0,tlist);
%     C(:,p)=OTOC(rho_0,H,V0,tlist);
%     C(:,p)=BiEETlist(rho_0,H,tlist,5);
end

% C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure(1)
plot(tlist/10,real(C)/real(C(1)))
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