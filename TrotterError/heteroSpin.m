% 20190613
% Compute the dynamics of two different spin species, one in dipolar state,
% the other in Z state. To see if they will reach thermal equilibrium.
% Model is a 1D 19F spin chain with one spin substitude by 31P
% Created by Pai Peng
N_atom=10;
N_P=5; % first to N_p spins are 31P
N_F=N_atom-N_P;
tic

bc='o';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
% model='Ising';
model=[1,1,-2]*8.18e-3;
% sym='kpz';
% sym='kp';
sym=[];
% tlist=(0:0.6281:10.5)/4;
tlist=logspace(-1,10,100);
% tlist=0:1:100;



Hz_P=OperatorClass(N_P,'z',1);
Hz_P.L=N_atom;
Hz_P.matrix={kron(Hz_P.matrix{1},speye(2^N_F))};

Hx_P=OperatorClass(N_P,'x',1);
Hx_P.L=N_atom;
Hx_P.matrix={kron(Hx_P.matrix{1},speye(2^N_F))};

Hy_P=OperatorClass(N_P,'y',1);
Hy_P.L=N_atom;
Hy_P.matrix={kron(Hy_P.matrix{1},speye(2^N_F))};

HXXYY_P=OperatorClass(N_P,[1,1,0],-1,bc,cplist);
HXXYY_P.L=N_atom;
HXXYY_P.matrix={kron(HXXYY_P.matrix{1},speye(2^N_F))};

Hz_F=OperatorClass(N_F,'z',1);
Hz_F.L=N_atom;
Hz_F.matrix={kron(speye(2^N_P),Hz_F.matrix{1})};

Hdipz_F=OperatorClass(N_F,model,-1,bc,cplist);
Hdipz_F.L=N_atom;
Hdipz_F.matrix={kron(speye(2^N_P),Hdipz_F.matrix{1})};

HXXYY_F=OperatorClass(N_F,[1,1,0],-1,bc,cplist);
HXXYY_F.L=N_atom;
HXXYY_F.matrix={kron(speye(2^N_P),HXXYY_F.matrix{1})};

HZZ=OperatorClass(N_atom,[0,0,-2],-1,bc,cplist);

Hff=LocalPauli(N_atom,N_P-1,'+')*LocalPauli(N_atom,N_P,'-')*LocalPauli(N_atom,N_P+1,'+')*LocalPauli(N_atom,N_P+2,'-');
Hff=(Hff+Hff')*1e-3;

rho_0=Hx_P+Hdipz_F;
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

H_int=OperatorClass(N_atom,model,-1,bc,cplist,2)+Hff;
% H_int=(HXXYY_F+HXXYY_P+HZZ)*8.18e-3;
H_int.symmetrize(sym);

Hz=2*pi*283*(Hz_F+17.2/40*Hz_P);
Hz.symmetrize(sym);
col=1;
    H=H_int+Hz;
%     H=U2H(H2U(0.25*H_int,3/180*pi)*H2U(0.5*Hz,3/180*pi),1);
    C=(abs(TPC(rho_0,H,Hx_P,tlist)).^2+abs(TPC(rho_0,H,Hy_P,tlist)).^2).^0.5;
%     C=TPC(rho_0,H,Hz_P,tlist);
%     C(:,p)=BiEETlist(rho_0,H,tlist,5);

% C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure(1)
semilogx(tlist,real(C))
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