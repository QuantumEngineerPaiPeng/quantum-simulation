% 20180805
% compute two-point correlator or OTOC involving an averaged operator 
% with trasverse field for different field strength
% and time
% Written using OperatorClass
% Created by Pai Peng
N_atom=10;

bc='p';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
model='dipolar';
sym='kpz';
% sym=[];
% tlist=(0:0.6281:10.5)/4;
tlist=logspace(1,5,50);
glist=linspace(0,3,7); % strength of the transverse field

rho_0=OperatorClass(N_atom,'z',1);
% rho_0=OperatorClass(N_atom);
% rho_0.matrix={zeros(2^N_atom,2^N_atom)};
% rho_0.matrix{1}(1,1)=1;
% V0=OperatorClass(N_atom,V0code,1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
V0=(copy(rho_0));

rho_0.symmetrize(sym);
V0.symmetrize(sym);

C=zeros(1,length(glist));

H_int=OperatorClass(N_atom,model,-1,bc,cplist,1);
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);
col=1;
parfor p=1:length(glist)
    g=glist(p)
    H=H_int-g*Hz;
    averho=AveOpe(rho_0, H, tlist);
%     C(:,p)=TPC(averho,H,V0,0);
    C(:,p)=OTOC(averho,H,V0,0);
end

C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure
hold on
plot(glist,C)
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
