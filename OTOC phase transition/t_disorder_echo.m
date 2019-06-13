%20180701
% compute two-point correlator or OTOC of an echo sequence with disorder 
% for different disorder strength and time
% Written using OperatorClass
% Created by Pai Peng
N_atom=8;
N_rand=20;
bc='p';
% cplist=[1,2,3].^-3;
cplist=1;
model='dipolar';
DisorderDir=[0,0,1]; %direction of the disorder field

tlist=logspace(-1,2,30);
wlist=linspace(0,5,6); % strength of the disorder field
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];

rho_0=OperatorClass(N_atom,'z',1);
% V0=OperatorClass(N_atom,'z',1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
V0=(copy(rho_0));

C=zeros(length(tlist),length(wlist),N_rand);

H_int=OperatorClass(N_atom,model,-1,bc,cplist,3);

for pp=1:N_rand
    Disorder=DisorderDir.'.*(rand(3,N_atom)*2-1);
    H_rand=randOpe(N_atom,Disorder);
    col=1;
    for w=wlist
        H1=H_int+w*H_rand;
        H2=H_int-w*H_rand;
        C(:,col,pp)=echoTPC(rho_0,H1,H2,V0,tlist);
        col=col+1;
    end
    pp
end

C=C/(N_atom*2^N_atom);
C_ave=real(sum(C,3))'/N_rand;

if N_rand==1
    figure(2)
    plot(tlist,C_ave)
else
    stdev=std(real(C),1,3)'/N_rand^0.5;
    figure
    hold on
    for p=1:length(wlist)
        errorbar(tlist,C_ave(p,:),stdev(p,:))
    end
end
