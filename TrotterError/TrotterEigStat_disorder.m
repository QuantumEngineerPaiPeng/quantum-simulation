% 20190704
% compute eigenvalue statistics for the trotter unitary. Use disorder to
% get rid of symmetries
% Created by Pai Peng
N_atom=13;
tic

bc='p';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
% model='Ising';
model='dipolar';
sym=[];
tlist=(40:5:65)/180*pi;

W=0.3;
N_rand=50;
DisorderDir=[1,0,0.3]; %direction of the disorder field
DisorderDir=DisorderDir/norm(DisorderDir);
seed=0;

R=zeros(length(tlist),N_rand);

H_int=OperatorClass(N_atom,model,-0.25,bc,cplist,1);

Hz=OperatorClass(N_atom,'z',1);
rng(seed)
randmat=rand(N_rand,N_atom);

parfor q=1:N_rand
    Disorder=DisorderDir.'*(randmat(q,:)*2-1);
    H_rand=randOpe(N_atom,Disorder);
    H1=H_int+W*H_rand;
    H1.diagonalize()
    rtemp=zeros(length(tlist),1);
    for p=1:length(tlist)
        t=tlist(p);
        
        U2=0.5*Hz;
%         U2.matrix={sparse(diag(exp(-1i*diag(U2.matrix{1})*t)))};
        U2.matrix={sparse(1:2^N_atom,1:2^N_atom,exp(-1i*diag(U2.matrix{1})*t)),2^N_atom,2^N_atom};
        U=H2U(H1,t)*U2;
        U.diagonalize()
        D=U.eigsys{1}.D;
        D=sort(real(log(D)/1i));
        delta=D(2:2^N_atom)-D(1:2^N_atom-1);
        for j=1:2^N_atom-2
            rtemp(p)=rtemp(p)+min([delta(j),delta(j+1)])/max([delta(j),delta(j+1)])/(2^N_atom-2);
        end
    end
    R(:,q)=rtemp;
end

R_ave=mean(R,2);

if N_rand==1
    figure
    plot(wlist,R_ave/(2*N_atom*2^N_atom))
else
    stdev=std(real(R),1,2)/N_rand^0.5;
    figure(3)
    errorbar(tlist,R_ave,stdev)
    hold on
end

% pfolder='~/Dropbox (MIT)/grad/research/codes/Tr/figure_data/';
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