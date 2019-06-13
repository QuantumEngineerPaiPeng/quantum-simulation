% 20180716
% function to get OTOC with averaged operator
function y=mblOTOCavg_func(B,N_phi,N_t1,N_t2)

Bb=reshape(real(B),N_phi,N_t1,N_t2,size(B,2));

for q=1:size(B,2)
    m=Bb(1,:,:,q);
    m=squeeze(m);
    n=sqrt(diag(m)*diag(m).');
    for p=1:N_phi
        m=Bb(p,:,:,q);
        m=squeeze(m);
        m=m./n;
        Bb(p,:,:,q)=m;
    end
end
Iq=fft(Bb,N_phi,1)/N_phi;
Iq=squeeze(mean(Iq,2));
Iq=squeeze(mean(Iq,2));
qsquare=(0:ceil(N_phi/2)).^2;
qsquare=[qsquare,(ceil(N_phi/2)-1:-1:1).^2];
y=qsquare*Iq;
end