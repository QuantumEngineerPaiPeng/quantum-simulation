% 20180629
% compute two-point correlator or OTOC with disorder for different disorder strength
% and time
% Written using OperatorClass
% Created by Pai Peng
tic
N_atom=10;
N_rand=10;
bc='p';
cplist=(1:5).^(-3);
% cplist=[1,2,3].^-3;
% cplist=1;
model=[0,1,-1];
DisorderDir=[1,0,0]; %direction of the disorder field
DisorderDir=DisorderDir/norm(DisorderDir);
seed=0;
savedata=0;

tlist=logspace(-2,2,50);
% wlist=1; % strength of the disorder field
% wlist=2*rand(10,1);
wlist=linspace(0,5,10);
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];
qu={'CL','Hamming',1};
% rho_0=(OperatorClass(N_atom,'x',1));
rho_0=(OperatorClass(N_atom,'y',1))-(OperatorClass(N_atom,'z',1));
% V0=OperatorClass(N_atom,'x',1);
% rho_0=LocalPauli(N_atom,3,'x');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=(OperatorClass(N_atom,'x',-1)+OperatorClass(N_atom,'y',1))*sqrt(1/2);
V0=(copy(rho_0));

C=zeros(length(tlist),length(wlist),N_rand);
% Ca=zeros(length(wlist),N_rand);
% Cweight=zeros(length(tlist),length(qu),length(wlist),N_rand);
% Ci=zeros(length(tlist),length(wlist),N_atom,N_rand);

H_int=OperatorClass(N_atom,model,-1,bc,cplist,3);
HB=(OperatorClass(N_atom,'z',1)+OperatorClass(N_atom,'y',1))*(1/sqrt(2));
rng(seed)
randmat=rand(N_rand,N_atom);
parfor pp=1:N_rand
    Disorder=DisorderDir.'*(randmat(pp,:)*2-1);
    H_rand=randOpe(N_atom,Disorder);

%     Citemp=zeros(length(tlist),length(wlist),N_atom);
    Ctemp=zeros(length(tlist),length(wlist));
%     Catemp=zeros(length(wlist),1);
%     Cweighttemp=zeros(length(tlist),3,length(wlist));
    for wp=1:length(wlist)
        w=wlist(wp);
        H=H_int+w*H_rand;
%         H.diagonalize()
%         averho=InfAveOpe(rho_0,H);
%         Catemp(wp)=TPC(averho,H,V0,0);
%         Ctemp(:,wp)=TPC(rho_0,H,V0,tlist);
        Ctemp(:,wp)=OTOC(rho_0,H,V0,tlist);
%         Cweighttemp(:,:,wp)=WeightsTlist(rho_0*sqrt(1/N_atom),H,tlist,qu);
%         for p=1:N_atom
%             rho0=LocalPauli(N_atom,p,'x');
%             Citemp(:,wp,p)=OTOC(rho0,H,V0,tlist);
%         end
    end
    C(:,:,pp)=Ctemp;
%     Ca(:,pp)=Catemp;
%     Ci(:,:,:,pp)=Citemp;
%     Cweight(:,:,:,pp)=Cweighttemp;
end

C_ave=real(sum(C,3))/N_rand/(2*N_atom*2^N_atom);
% Caave=mean(Ca,2)/2^N_atom/N_atom;
% Cweightave=mean(Cweight,4);

% sCi=squeeze(sum(Ci,3));
% if length(wlist)==1
%     sCiave=mean(sCi,2)/(2*N_atom*2^N_atom);
% else
%     sCiave=mean(sCi,3)/(2*N_atom*2^N_atom);
% end
if N_rand==1
    figure(2)
    plot(tlist,C_ave.')
else
    Cstd=std(real(C),1,3)/N_rand^0.5/(2*N_atom*2^N_atom);
%     Castd=std(Ca,1,2)/N_rand^0.5/2^N_atom/N_atom;
%     Cweightstd=std(Cweight,1,4)/(N_rand)^0.5;
%     if length(wlist)==1
%         sCistd=std(sCi,1,2)/(N_rand)^0.5/(2*N_atom*2^N_atom);
%     else
%         sCistd=std(sCi,1,3)/(N_rand)^0.5/(2*N_atom*2^N_atom);
%     end
    
    figure
    hold on
%     errorbar(wlist,Caave,Castd)
%     set(gca,'xscale','log')
    for p=1:length(wlist)
        errorbar(tlist,C_ave(:,p),Cstd(:,p))
        
%         errorbar(tlist,sCiave(:,p),sCistd(:,p))
%         for q=1:length(qu)
%             errorbar(tlist,Cweightave(:,q,p),Cweightstd(:,q,p))
%         end
        set(gca,'xscale','log')
    end

    xlabel('t')
    ppStyle(30,2,10)
end

if savedata
    if isunix
        pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
    else
        pfolder='C:\\Users\\Pai\\Dropbox (MIT)\\grad\\research\\codes\\OTOC phase transition\\figure_data\\';
    end
    foldername=mfilename;
    mkdir(pfolder,foldername)
    if isunix
        filetosave=[pfolder,foldername,'/'];
    else
        filetosave=[pfolder,foldername,'\\'];
    end
    for p = fix(clock)
        filetosave=[filetosave,num2str(p),'_'];
    end
    filetosave=[filetosave,'.mat'];
    save(filetosave)
    varname={'file','N_atom','N_rand','model','wlist','bc','cplist','DisorderDir','seed'};
    varvalue={filetosave,N_atom,N_rand,model,wlist,bc,cplist,DisorderDir,seed};
    writelog(pfolder,varname,varvalue)
end

toc