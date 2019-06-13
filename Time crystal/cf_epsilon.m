%20170913
%crystalline fraction v.s. epsilon in various time window
N_atom=9;
N_rand=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]
tau=0;

N=2^6;
n_wins=8;
rho_0=Sigma_y;

min=0;
max=0.15;
n_step=9;
step=(max-min)/n_step;

x=min:step:max;

cf=zeros(n_wins,n_step+1,N_rand);
F=zeros(n_wins,N/n_wins);

H_int=zeros(2^N_atom);

for p=1:N_atom-1 %nearest neighbor
    H_int=H_int+(-1)*...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
end
% H_int=H_int+(-1)*... %periodic boundary condition
%     (Sigma_i{2}{N_atom}*Sigma_i{2}{1}-...
%     (Sigma_i{3}{N_atom}*Sigma_i{3}{1}+Sigma_i{1}{N_atom}*Sigma_i{1}{1})/2);
% for p=1:N_atom-2 %next-nearest neighbor
%     H_int=H_int+(-1/8)*...
%         (Sigma_i{2}{p}*Sigma_i{2}{p+2}-...
%         (Sigma_i{3}{p}*Sigma_i{3}{p+2}+Sigma_i{1}{p}*Sigma_i{1}{p+2})/2);
% end
H_rand=zeros(2^N_atom);

for pp=1:N_rand
    count=1;
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(rand-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(rand-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(rand-0.5)*Sigma_i{3}{p};
    end
    
    H=H_int+H_rand;
    
    U_H=exp(1)^(-1i*H*tau);
    
    for epsilon=min:step:max
        U_p=exp(1)^(-1i*Sigma_x*pi/2*(1+epsilon));
        
        U=U_p*U_H;
        
        %     t=0:N-1;
        %     t=t*tau;
        rhot=cell(1,N);
        rhot{1}=rho_0;
        M_y=zeros(1,N);
        %     M_x=zeros(1,N);
        
        for p=1:N
            rhot{p+1}=U*rhot{p}*U^(-1);
            M_y(p)=trace(rhot{p+1}*Sigma_y/2^N_atom);
            %         M_x(p)=trace(rhot{p+1}*Sigma_x/2^N_atom);
        end
        
        for p=1:n_wins
            F(p,:)=abs(fft(real(M_y(1+(p-1)*N/n_wins:p*N/n_wins))));
        end
        
        %cf(:,count,pp)=F(:,N/2/n_wins+1)./sum(F,2);
        cf(:,count,pp)=F(:,N/2/n_wins+1).^2./sum(F.^2,2);
        count=count+1;
    end
    pp
%     figure(1)
%     hold on
%     plot(x,cf)
end
figure(1)
y=sum(cf,3)/N_rand;
stdev=std(cf,1,3)/N_rand^0.5;
hold on
for p=1:n_wins
    errorbar(x,y(p,:),stdev(p,:))
end