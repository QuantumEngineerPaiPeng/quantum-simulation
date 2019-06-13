%20170912
%evolve an operator
N_atom=8;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

J=1;
g=1;

% tau=0.5;
% epsilon=0.1/pi;
% N=200;
% tlist=logspace(1,10,30);
tlist=0:0.1:3;
% rho_0=diag(zeros(1,2^N_atom));
% rho_0(1,1)=1;
rho_0=Sigma_z;
cplist=[1,1/8,1/27];
bc='p';

% frame=1:(N+1);
% frame=mod(frame,2)*2-1;
H_int=Hamiltonian(N_atom,bc,1,'dipolar',1);

H=H_int+g*Sigma_z;

% U_H=expm(-1i*H*tau);
% U_p=expm(-1i*pi*(1+epsilon)*Sigma_z/2);
% U=U_H*U_p;

% t=0:N;
% t=t*tau;
% rhot=rho_0;
% M_y=zeros(1,N+1);
M_x=zeros(1,length(tlist));
% M_z=zeros(1,N+1);
% M_x(1)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
[V,D]=eig(H);
count=1;
for t=tlist
    U=V*diag(exp(1i*diag(D)*t))*V';
    rhot=U*rho_0*U';
%     M_y(p)=trace(rhot{p}*Sigma_y/(N_atom*2^N_atom));
    M_x(count)=trace(rhot*rho_0);%/((N_atom-1)*2^N_atom));
%     M_z(p)=trace(rhot{p}*Sigma_z/(N_atom*2^N_atom));
count=count+1;
    
end

% F=abs(fft(real(M_y)));

figure(1)
hold on
% plot(tlist,M_x,'LineWidth',2)
plot(real(M_x))

% figure(2)
% plot((0:N-1)/N/tau,F)