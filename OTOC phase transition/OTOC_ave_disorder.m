%20180119
%compute OTOC v.s w with one of the operators being long-time averaged
N_atom=9;
N_rand=50;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[1,1,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

wlist=logspace(-1,2,20);
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];

rho_0=(Sigma_x+Sigma_z)/sqrt(2);
% rho_0=Sigma_x;
V0=rho_0;

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
C=zeros(length(wlist),N_rand);

for p=1:N_atom-1
%     H_int=H_int+(-1)*...
%         (Sigma_i{1}{p}*Sigma_i{1}{p+1}-...
%         (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{2}{p}*Sigma_i{2}{p+1})/2);
    %     H_int=H_int+(-1)*...
    %         (Sigma_i{2}{p}*Sigma_i{2}{p+1});
    %     H_int=H_int+...
    %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        H_int=H_int+(-1)*...
        (Sigma_i{1}{p}*Sigma_i{1}{p+1}-Sigma_i{3}{p}*Sigma_i{3}{p+1});
end
% H_int=H_int+(-1)*...
%     (Sigma_i{1}{N_atom}*Sigma_i{1}{1}-...
%     (Sigma_i{3}{N_atom}*Sigma_i{3}{1}+Sigma_i{2}{N_atom}*Sigma_i{2}{1})/2);
rng(ss)
for pp=1:N_rand
    Disorder=rand(N_atom,1)*2-1;
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(p)*(Sigma_i{1}{p}+Sigma_i{3}{p})/sqrt(2);
    end
    col=1;
    for w=wlist
        H=H_int+w*H_rand;
        [V,D]=eig(H);
        row=1;
        rho_ave=V'*rho_0*V;
        C(col,pp)=-trace((diag(diag(rho_ave))*rho_ave-rho_ave*diag(diag(rho_ave)))^2);
        col=col+1;
    end
    pp
end
C_ave=real(sum(C,2))'/N_rand;

if N_rand==1
    figure(2)
    plot(tlist,C_ave/(2*N_atom*2^N_atom))
else
    stdev=std(real(C),1,2)'/N_rand^0.5;
    hold on
    figure(1)
    errorbar(wlist.^-1,C_ave,stdev)
end
