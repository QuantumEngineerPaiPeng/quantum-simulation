% 20190607
% compute fidelity of trotterized Hdipz+X
% Written using OperatorClass
% Created by Pai Peng
N_atom=12;
tic

bc='p';
cplist=[1,2,3].^-3;
model=[-0.5,-0.5,1];
% sym='kpz';
sym='kp';
% sym=[];
% tlist=(0:0.6281:10.5)/4;
% tlist=logspace(-1,5,50);
tlist=0:5:500;

% philist=30:30:150;
philist=26:2:34;

rho_0=OperatorClass(N_atom,'y',1);
rho_0.symmetrize(sym);

V0=rho_0;

C=zeros(length(tlist),length(philist));

H_int=OperatorClass(N_atom,model,-1,bc,cplist,2);
H_int.symmetrize(sym);

Hx=OperatorClass(N_atom,'x',1);
Hx.symmetrize(sym);
Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);
col=1;

Hb=-0.25*H_int-0.5*Hx;
Hb.diagonalize()

parfor pp=1:length(philist)
    phi=philist(pp);
    Uf=H2U(0.25*H_int,30/180*pi)*H2U(0.5*Hx,phi/180*pi);
    Hf=U2H(Uf,30/180*pi);
    
    ctemp=zeros(length(tlist),1);
    for p=1:length(tlist)
        t=tlist(p);
        Uft=H2U(Hf,t);
        Ubt=H2U(Hb,t);
        ctemp(p)=overlap(Ubt*Uft*rho_0*Uft'*Ubt',V0);
        %     C(:,p)=OTOC(rho_0,H,V0,tlist);
        %     C(:,p)=BiEETlist(rho_0,H,tlist,5);
    end
    C(:,pp)=ctemp;
end

% C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure
plot(tlist,real(C)/(N_atom*2^N_atom))
xlabel('Jt')
ylabel('Fidelity')
% ppStyle(30,2,10)
lgd={};
for p=1:length(philist)
    lgd{end+1}=num2str(philist(p));
end
legend(lgd)
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