%20170912
%simulate time crystal dynamics [pi exp(iHt)]^n
N_atom=8;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

J=1;
g=0.;

tau=0.5;
epsilon=0.1/pi;
N=10;
rho_0=Sigma_x;
cplist=[1,1/8,1/27];
bc='p';

frame=1:(N+1);
frame=mod(frame,2)*2-1;
H_int=-Hamiltonian(N_atom,bc,cplist,'dipolar',1);

H=H_int+g*Sigma_z;

U_H=expm(-1i*H*tau);
U_p=expm(-1i*pi*(epsilon)*Sigma_z/2);
U=U_p*U_H;

t=0:N;
t=t*tau;
rhot=rho_0;
% M_y=zeros(1,N+1);
M_x=zeros(1,N+1);
% M_z=zeros(1,N+1);
M_x(1)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
for p=2:N+1
    rhot=U*rhot*U';
%     M_y(p)=trace(rhot{p}*Sigma_y/(N_atom*2^N_atom));
    M_x(p)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
%     M_z(p)=trace(rhot{p}*Sigma_z/(N_atom*2^N_atom));
    
end

% F=abs(fft(real(M_y)));

figure(2)
hold on
plot(0:N,M_x,'LineWidth',2)
% figure(2)
% plot((0:N-1)/N/tau,F)