% 20180803
% compute two-point correlator or OTOC involving an averaged operator 
% with disorder for different disorder strength and time
% Written using OperatorClass
% Created by Pai Peng
N_atom=8;
N_rand=50;
bc='o';
% cplist=[1,2,3].^-3;
cplist=1;
model='dipolar';
DisorderDir=[0,0,1]; %direction of the disorder field
DisorderDir=DisorderDir/norm(DisorderDir);

tlist=4.2:0.2:5;
wlist=logspace(0,2,30); % strength of the disorder field
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];

rho_0=OperatorClass(N_atom,'z',1);
V0=OperatorClass(N_atom,'z',1);
% rho_0=LocalPauli(N_atom,1,'x');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
% V0=(copy(rho_0));
% V0=

C=zeros(length(wlist),N_rand);

H_int=OperatorClass(N_atom,model,-1,bc,cplist,3);

for pp=1:N_rand
    Disorder=DisorderDir.'*(rand(1,N_atom)*2-1);
    H_rand=randOpe(N_atom,Disorder);
%     H_rand=OperatorClass(N_atom,'x',1);

    parfor wp=1:length(wlist)
        w=wlist(wp);
        H=H_int+w*H_rand;
        H.diagonalize()
%         averho=AveOpe(rho_0, H, tlist);
        averho=InfAveOpe(rho_0,H);
        C(wp,pp)=TPC(averho,H,V0,0);
%         C(wp,pp)=TwoNorm(averho);
%         C(wp,pp)=OTOC(averho,H,V0,0);
    end
    pp
end
C=C/(2*N_atom*2^N_atom);
C_ave=real(sum(C,2))'/N_rand;


if N_rand==1
    figure(3)
    plot(wlist,C_ave)
else
    stdev=std(real(C),1,2)'/N_rand^0.5;
    figure(3)
    hold on
    errorbar(wlist,C_ave,stdev)
    xlabel('w')
    ppStyle(30,2,10)
end
