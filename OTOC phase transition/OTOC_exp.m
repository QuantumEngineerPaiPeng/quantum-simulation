%20180222
%anlyse the otoc data from experiments
%data is read by data_analysis file to variable 'B'
N_phi=6;% number of encoding angles
N_cycle=4;% number of phase cycles
N_t=12;
Bb=reshape(B,N_phi,N_t,N_cycle,size(B,2));
ERR=reshape(err,N_phi*N_cycle*N_t*size(B,2),1);
% VAR=sparse(diag(ERR.^2));
Bb=real((mean(Bb,3)));
Bb=Bb./Bb(1,:,:,:);
Bb=squeeze(Bb);
Iq=fft(Bb,N_phi,1)/N_phi;
Iq=abs(Iq);
Iq=Iq./sum(Iq,1);
qsquare=(0:ceil(N_phi/2)).^2;
qsquare=[qsquare,(ceil(N_phi/2)-1:-1:1).^2];
OTOC=qsquare*(reshape(Iq,N_phi,N_t*(size(B,2))));
OTOC=reshape(OTOC,N_t,(size(B,2)));
% figure(1)
% hold on
% plot(1:size(Bb,3),OTOC(3:3:end,:),'r','LineWidth',2)
% plot(1:size(Bb,3),OTOC,'LineWidth',2)