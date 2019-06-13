%20170913
%time evolution averaged over disorder realizations
N_atom=8;
N_rand=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]
tau=0.2945;

N=2^8;

epsilon=0.1/pi;

rho_0=Sigma_z;
% rho_0=zeros(2^N_atom);
% rho_0(1,1)=1;
% rho_0=exp(1)^(1i*Sigma_x*pi/4)*rho_0*exp(1)^(-1i*Sigma_x*pi/4);
U_p=exp(1)^(-1i*Sigma_x*pi/2*(1+epsilon));

%cf=zeros(n_wins,n_step+1,N_rand);
%F=zeros(n_wins,N/n_wins);
M_y=zeros(N_rand,N);
RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);

for p=1:N_atom-1
    H_int=H_int+(-1)*...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}-...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
    %     H_int=H_int+(-2)*...
    %         (Sigma_i{2}{p}*Sigma_i{2}{p+1});
end

H_rand=zeros(2^N_atom);

for pp=1:N_rand
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(RandomM(pp,1,p)-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(RandomM(pp,2,p)-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(RandomM(pp,3,p)-0.5)*Sigma_i{3}{p};
    end
    
    H=H_int+H_rand;
    
    U_H=exp(1)^(-1i*H*tau);
    
    U=U_p*U_H;
    
    %     t=0:N-1;
    %     t=t*tau;
    rhot=cell(1,N);
    rhot{1}=rho_0;
    %M_y=zeros(1,N);
    %     M_x=zeros(1,N);
    
    for p=1:N
        rhot{p+1}=U*rhot{p}*U^(-1);
        M_y(pp,p)=trace(rhot{p+1}*Sigma_z/(N_atom*2^N_atom))*(-1)^(p);
        %         M_x(p)=trace(rhot{p+1}*Sigma_x/2^N_atom);
    end
    
    %     for p=1:n_wins
    %         F(p,:)=abs(fft(real(M_y(1+(p-1)*N/n_wins:p*N/n_wins))));
    %     end
    
    pp
    %     figure(1)
    %     hold on
    %     plot(x,cf)
end
Y=sum(real(M_y),1)/N_rand;
stdev=std(real(M_y),1,1)/N_rand^0.5;
if N_rand==1
    figure(1)
    plot(1:N,Y)
else
    figure(1)
    hold on
    color='r';
    plot(1:N,Y,color)
    plot(1:N,Y+stdev,color)
    plot(1:N,Y-stdev,color)
    hold off
end