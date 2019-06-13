%20180710
%anlyse the otoc data from experiments
%data is read by data_analysis file to variable 'B'
% Error analysis is included
N_phi=8;% number of encoding angles
N_cycle=1;% number of phase cycles
N_t=15;
Bb=reshape(real(B),N_phi,N_t,N_cycle,size(B,2));
ERR=reshape(err,N_phi,N_t,N_cycle,size(B,2),1);
otoc=zeros(N_t,size(B,2));
ebar=zeros(N_t,size(B,2));
for t=1:N_t
    for p=1:size(B,2)
                bb=squeeze(Bb(:,t,:,p));
                bb=real(mean(bb,2));

                err1=squeeze(ERR(:,t,:,p));
                err1=reshape(err1,N_phi*N_cycle,1);

        VAR=diag(err1.^2);

                VAR=reshape(VAR,N_phi,N_cycle,N_phi,N_cycle);
                VAR=permute(VAR,[2,1,3,4]);
                T=ones(1,N_cycle)/N_cycle;
                VAR=reshape(VAR,N_cycle,N_phi*N_phi*N_cycle);
                VAR=T*VAR;
                VAR=reshape(VAR,N_phi*N_phi,N_cycle);
                VAR=VAR*T.';
                VAR=VAR/(bb(1)^2);
                bb=bb/bb(1);
        VAR=reshape(VAR,N_phi,N_phi);
        T=exp(-2*pi*1i/N_phi*((0:N_phi-1).'*(0:N_phi-1)))/N_phi;
        iq=T*bb;
        T1=diag(real(iq)./(abs(iq)))*real(T)+diag(imag(iq)./(abs(iq)))*imag(T);
        iq=abs(iq);
        T2=diag(sum(iq)^(-1)*ones(N_phi,1))-kron(sum(iq)^(-2)*iq,ones(1,N_phi));
        iq=iq/sum(iq);
        qsquare=(0:ceil(N_phi/2)).^2;
        qsquare=[qsquare,(ceil(N_phi/2)-1:-1:1).^2];
        otoc(t,p)=qsquare*iq;
        Tt=qsquare*T2*T1;
        ebar(t,p)=sqrt(Tt*VAR*Tt.');
    end
end

% figure(4)
% hold on
% plot(1:size(B,2),otoc(3:3:end,:),'b','LineWidth',2)
% [6,9,12,15,18,21,24,30,36,42,48]
% for p=1:size(otoc,1)
% errorbar([6,12,18,24,30,36,42,48],otoc(p,:),ebar(p,:),':')
% end