%20180222
%anlyse the otoc data from experiments
%data is read by data_analysis file to variable 'B'
N_phi=6;% number of encoding angles
N_cycle=4;% number of phase cycles
N_t=12;
Bb=reshape(real(B),N_phi,N_t,N_cycle,size(B,2));
ERR=reshape(err,N_phi*N_cycle*N_t*size(B,2),1);
VAR=sparse(diag(ERR.^2));
% Bb=real((mean(Bb,3)));
Bb=permute(Bb,[3,1,2,4]);
Bb=reshape(Bb,N_cycle,N_phi*N_t*size(B,2));
VAR=reshape(VAR,N_cycle,N_phi*N_t*size(B,2)*N_cycle*N_phi*N_t*size(B,2));
T=ones(1,N_cycle)/N_cycle;
Bb=T*Bb;
VAR=T*VAR;
VAR=reshape(VAR,N_phi*N_t*size(B,2),N_cycle,N_phi*N_t*size(B,2));
VAR=permute(VAR,[1,3,2]);
VAR=reshape(VAR,N_phi*N_t*size(B,2)*N_phi*N_t*size(B,2),N_cycle);
VAR=VAR*T';
% Bb=Bb./Bb(1,:,:,:);
Bb=reshape(Bb,N_phi,N_t*size(B,2));
T=kron(ones(N_phi,1),Bb(1,:).^-1);
T=reshape(T,N_phi*N_t*size(B,2),1);
T=sparse(diag(T));
Bb=reshape(Bb,N_phi*N_t*size(B,2),1);
Bb=T*Bb;
VAR=reshape(VAR,N_phi*N_t*size(B,2),N_phi*N_t*size(B,2));
VAR=T*VAR*T';
Bb=reshape(Bb,N_phi,N_t*size(B,2));
T=exp(-2*pi*1i/N_phi*((0:N_phi-1).'*(0:N_phi-1)))/N_phi;
% Iq=fft(Bb,N_phi,1)/N_phi;
Iq=T*Bb;
T=kron(T,speye(N_t*size(B,2)));
vIq=reshape(Iq,N_phi*N_t*size(B,2),1);
T=diag(real(vIq)./(abs(vIq)))*real(T)+diag(imag(vIq)./(abs(vIq)))*real(T);
qsquare=(0:ceil(N_phi/2)).^2;
qsquare=[qsquare,(ceil(N_phi/2)-1:-1:1).^2];
Iq=abs(Iq);
A=sum(Iq,1);
Iq=Iq./sum(Iq,1);
T1=kron(ones(N_phi,1),sum(Iq,1));
T1=reshape(T1,N_phi*N_t*size(B,2),1).^(-1);
T11=diag(T1);

otoc=qsquare*Iq;
otoc=reshape(otoc,N_t,(size(B,2)));
VAR=reshape(VAR,N_phi*N_t*size(B,2),N_phi*N_t*size(B,2));
VAR=T*VAR*T.';
% VAR=qsquare*T*VAR;
% VAR=T*VAR;
% VAR=reshape(VAR,N_t*size(B,2),N_phi,N_t*size(B,2));
VAR=reshape(VAR,N_phi*N_t*size(B,2),N_phi,N_t*size(B,2));
VAR=permute(VAR,[1,3,2]);
% VAR=reshape(VAR,N_t*size(B,2)*N_t*size(B,2),N_phi);
VAR=reshape(VAR,N_phi*N_t*size(B,2)*N_t*size(B,2),N_phi);
% VAR=VAR*T'*qsquare';
VAR=VAR*T';
VAR=reshape(VAR,N_t*size(B,2),N_t*size(B,2));
STD=real(reshape(diag(VAR).^0.5,N_t,size(B,2)));

figure(1)
hold on
% plot(1:size(B,2),otoc(3:3:end,:),'b','LineWidth',2)
plot(1:size(B,2),otoc,'LineWidth',2)