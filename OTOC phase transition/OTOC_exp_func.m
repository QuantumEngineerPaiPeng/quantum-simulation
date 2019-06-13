%20180716
% function to anlyse the otoc data from experiments
%data is read by data_analysis file to variable 'B'
% function [y,Iq]=OTOC_exp_func(B,N_phi,N_cycle,N_t)
N_phi=8;% number of encoding angles
N_cycle=1;% number of phase cycles
N_t=8;
if length(size(B))==3
    B=reshape(B,size(B,1)*size(B,2),size(B,3)); % for MBL OTOC
end
Bb=reshape(B,N_phi,N_t,N_cycle,size(B,2));
Bb=real((mean(Bb,3)));
fidelity=squeeze(Bb(1,:,:));
Bb=Bb./Bb(1,:,:,:);
Bb=squeeze(Bb);
Iq=fft(Bb,N_phi,1)/N_phi;
Iq=abs(Iq);
Iq=Iq./sum(Iq,1);
qsquare=(0:ceil(N_phi/2)).^2;
qsquare=[qsquare,(ceil(N_phi/2)-1:-1:1).^2];
otoc=qsquare*(reshape(Iq,N_phi,N_t*(size(B,2))));
y=reshape(otoc,N_t,(size(B,2)));
% end