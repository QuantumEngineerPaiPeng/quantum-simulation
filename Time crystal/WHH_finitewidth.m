%20180322
%simulate WAHUHA sequence with finite pulse width
N_atom=8;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

J=8.18e-3;

cycle=30;
epsilon=0.02;
N=100;
rho_0=Sigma_x;
bc='p';
% cplist=[1,1/8,1/27];
cplist=1;

H_int=zeros(2^N_atom);
for pp=1:length(cplist)
    for p=1:N_atom-1
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+J*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
        %     H_int=H_int+(-1)*...
        %         (Sigma_i{1}{p}*Sigma_i{1}{p+1});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end

px=expm(1i*(H_int-pi/(4*p1)*Sigma_x)*p1);% here are U^\dag
py=expm(1i*(H_int-pi/(4*p1)*Sigma_y)*p1);
pxb=expm(1i*(H_int+pi/(4*p1)*Sigma_x)*p1);
pyb=expm(1i*(H_int+pi/(4*p1)*Sigma_y)*p1);

tau=cycle/6;

t1=tau-p1/2;
t2=tau-p1;

Uf=V*diag(exp(diag(1i*D*t1)))*V'*px*V*diag(exp(diag(1i*D*t2)))*V'*py*...
    V*diag(exp(diag(1i*2*D*t1)))*V'*pyb*V*diag(exp(diag(1i*D*t2)))*V'*...
    pxb*V*diag(exp(diag(1i*D*t1)))*V';

t=0:N;
t=t*cycle;
rhot=rho_0;
% M_y=zeros(1,N+1);
M_x=zeros(1,N+1);
% M_z=zeros(1,N+1);
M_x(1)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
for p=2:N+1
    rhot=Uf*rhot*Uf';
%     M_y(p)=trace(rhot{p}*Sigma_y/(N_atom*2^N_atom));
    M_x(p)=trace(rhot*Sigma_x/(N_atom*2^N_atom));
%     M_z(p)=trace(rhot{p}*Sigma_z/(N_atom*2^N_atom));
end

% F=abs(fft(real(M_y)));

figure(1)
hold on
plot(0:N,M_x,'LineWidth',2)
% figure(2)
% plot((0:N-1)/N/tau,F)