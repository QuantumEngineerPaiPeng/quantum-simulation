% 20180724
% for an initial state rho_0, evaluate the corresponding canonical ensemble
% averages, assuming ETH
N_atom=14;

bc='p';
% bc='o';
cplist=[1,2,3].^-3;
% cplist=1;
% cplist=(1:N_atom-1).^-3;
model='dipolar';
sym=['kpz'];
% sym=[];
% tlist=(0:0.6281:10.5)/4;
% tlist=logspace(-1,3,30);
glist=1; % strength of the transverse field
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];
% rho_0code='z';
V0code='z';
% rho_0=OperatorClass(N_atom,rho_0code,1);
rho_0=OperatorClass(N_atom);
rho_0.matrix={zeros(2^N_atom,2^N_atom)};
rho_0.matrix{1}(1,1)=1;
V0=OperatorClass(N_atom,V0code,1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
% V0=(copy(rho_0));

rho_0.symmetrize(sym);
V0.symmetrize(sym);

Ob=zeros(1,length(glist));

H_int=OperatorClass(N_atom,model,-1,bc,cplist,1);
H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);
col=1;
for p=1:length(glist)
    g=glist(p);
    H=1/4*H_int-g/2*Hz;
    E=overlap(H,rho_0);
    fun=@(x) CanEnsEnergy(x,H)-E;
    beta=fsolve(fun,0);
    Z=CanEnsPartition(beta,H);
    rho=1/Z*H2U(-1i*beta*H,1);
    Ob(p)=overlap(rho,V0);
end

% C=C/(2*N_atom*2^N_atom);
% C=C/N_atom;
% C_ave=mean(C,1);
% 
% stdev=std(C,1)/sqrt(length(tlist));
figure
hold on
plot(glist,Ob)
% xlabel('g')
ppStyle(30,2,10)