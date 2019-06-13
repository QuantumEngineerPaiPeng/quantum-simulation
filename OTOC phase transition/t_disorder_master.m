% 20180803
% compute two-point correlator or OTOC with disorder for different disorder strength
% and time
% Consider the full master sequence
% Written using OperatorClass
% Created by Pai Peng
N_atom=8;
N_rand=50;
bc='p';
cplist=[1,2,3].^-3;
% cplist=1;
model='dipolar';

eps=6;
eps=eps/180*pi;

cycle=100;
tau=cycle/24;
J=8.18e-3;

a=0;
b=0.6;
c=0;
u=-0.2;
v=0.2;
w=0;

p1=1;

tlist=1:64;
wlist=linspace(0,2,5); % strength of the disorder field
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];

rho_0=OperatorClass(N_atom,'x',1);
V0=OperatorClass(N_atom,'x',1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
% V0=(copy(rho_0));
% V0=

C=zeros(length(tlist),length(wlist),N_rand);

H_int=OperatorClass(N_atom,model,-J,bc,cplist,1);
Z=OperatorClass(N_atom,'z',1);

phase=[0,90,90,0,0,90,90,0,180,270,270,180,180,270,270,180];
angle=pi/2;

t1=tau*(1-v+w+c)-p1/2;
t11=tau*(1-v+w-c)-p1/2;
t2=tau*(1-u+v+b)-p1;
t22=tau*(1-u+v-b)-p1;
t3=tau*(1+u-w+c)-p1/2;
t33=tau*(1+u-w-c)-p1/2;

delay=[t1,t2,2*t3,t22,2*t11,t2,2*t33,t22,2*t1,t22,2*t33,t2,2*t11,t22,2*t3,t2,t1];

for pp=1:N_rand
    Disorder=[0,0,1].'*(rand(1,N_atom)*2-1);
    H_rand=randOpe(N_atom,Disorder);
	H=H_int+H_rand;
    col=1;
    for w=wlist
        Uf=UFloquet(phase,angle,delay,H,p1)*H2U(Z,eps/2);
        Hf=U2H(Uf,1);
%         H=H_int+w*H_rand;
%         H.diagonalize()
        C(:,col,pp)=TPC(rho_0,Hf,V0,tlist);
%         C(:,col,pp)=OTOC(rho_0,H,V0,tlist);
        col=col+1;
    end
    pp
end
C=C/(2*N_atom*2^N_atom);
C_ave=real(sum(C,3))'/N_rand;


if N_rand==1
    figure(3)
    plot(tlist,C_ave)
else
    stdev=std(real(C),1,3)'/N_rand^0.5;
    figure(3)
    hold on
    for p=1:length(wlist)
        errorbar(tlist,C_ave(p,:),stdev(p,:))
    end
    xlabel('t')
    ppStyle(30,2,10)
end
