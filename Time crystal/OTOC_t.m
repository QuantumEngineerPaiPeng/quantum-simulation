%20180102
%compute OTOC v.s t
N_atom=8;
N_rand=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

rho_0=Sigma_z;

tlist=0:0.1:1;
N_t=length(tlist);
g=0.5;

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
C=zeros(N_t,N_rand);

for p=1:N_atom-1
%     H_int=H_int+(-1)*...
%         (Sigma_i{1}{p}*Sigma_i{1}{p+1}-...
%         (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{2}{p}*Sigma_i{2}{p+1})/2);
            H_int=H_int+(-1)*...
                (Sigma_i{1}{p}*Sigma_i{1}{p+1});
%     H_int=H_int+...
%         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
    
end

for pp=1:N_rand
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(rand-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(rand-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(rand-0.5)*Sigma_i{3}{p};
    end
    
    count_t=1;
    for t=tlist
        H=H_int+H_rand+g*Sigma_z;
        U=exp(1)^(-1i*H*t);
        rhot=U*rho_0*(inv(U));
        C(count_t,pp)=-trace((rhot*Sigma_x-Sigma_x*rhot)^2);
        count_t=count_t+1;
    end
    
end
C_ave=real(sum(C,2))'/N_rand;

if N_rand==1
    figure(1)
    plot(tlist,C_ave)
else
    stdev=std(real(C),1,3)'/N_rand^0.5;
    figure(1)
    hold on
    for p=1:N_t
        errorbar(tlist,C_ave(:,p),stdev(:,p))
    end
end
