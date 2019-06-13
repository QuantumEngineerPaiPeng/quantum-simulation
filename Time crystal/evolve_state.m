%20170912
%evolve a pure state
N_atom=8;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

J=1;
g=16;

% tau=0.5;
% epsilon=0.1/pi;
% N=200;
tlist=logspace(-2,1,50);
% tlist=linspace(0,10,20);
psi0=zeros(2^N_atom,1);
psi0(1)=1;
cplist=[1,1/8,1/27];
bc='p';

% frame=1:(N+1);
% frame=mod(frame,2)*2-1;
H_int=-Hamiltonian(N_atom,bc,1,'dipolar',1);

H=H_int+g*Sigma_z;

% U_H=expm(-1i*H*tau);
% U_p=expm(-1i*pi*(1+epsilon)*Sigma_z/2);
% U=U_H*U_p;

% t=0:N;
% t=t*tau;
% rhot=rho_0;
% M_y=zeros(1,N+1);
y=zeros(1,length(tlist));
% M_z=zeros(1,N+1);
% M_x(1)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
[V,D]=eig(H);
count=1;
for t=tlist
    U=V*diag(exp(1i*diag(D)*t))*V';
    psi=U*psi0;
    Gs=reshape(psi,[2^(N_atom/2),2^(N_atom/2)]);
    s=svd(Gs);
    y(count)=-2*log(s')*(s.^2);%/((N_atom-1)*2^N_atom));
    
count=count+1;
    
end

% F=abs(fft(real(M_y)));

figure(2)

semilogx(tlist,abs(y).^2,'LineWidth',2)
hold on
% figure(2)
% plot((0:N-1)/N/tau,F)