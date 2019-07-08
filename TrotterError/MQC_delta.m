% 20190701
% Calculate the MQC of dipolar Z state
% created by Pai Peng
N_atom=10;
tic

bc='p';
% cplist=[1,2,3].^-3;
cplist=1;
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
Hy=OperatorClass(N_atom,'y',1);
Hy.symmetrize(sym);
U45=H2U(Hx,pi/4/2);

V0=OperatorClass(N_atom,'x',1);
% V0=OperatorClass(N_atom,[0,1,-1],-1,bc,cplist);
V0.symmetrize(sym);

C=zeros(length(tlist),q_max);

H_int=OperatorClass(N_atom,model,-1,bc,cplist);
H_int.symmetrize(sym);
H_int.diagonalize()

Ufid=H2U(H_int,0.3); % unitary after the final pi/4

V0=U45'*Ufid'*V0*Ufid*U45;

for q=1:q_max
    Uphi=H2U(Hx,(q-1)*2*pi/q_max/2);
    C(:,q)=TPC(Uphi*rho_0*Uphi',H_int,V0,tlist);
end

for p=1:length(tlist)
    C(p,:)=C(p,:)/sqrt(C(p,:)*C(p,:)');
end
% C=C/(2*N_atom*2^N_atom);

% figure
% plot(((1:q_max)-1)/q_max,real(C))
% hold on
% xlabel('\phi (2\pi)')
% ppStyle(30,2,10)

Iq=fft(C,[],2);
figure(4)
plot(tlist/2,abs(Iq(:,[1,3]))*2)
xlabel('\delta')
ylabel('I_q')
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