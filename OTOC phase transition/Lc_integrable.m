%20180121
%Compute OTOC for intergrable model
N_spin=10;
J=1;
g=0.5;
%t=0.16*12;
tlist=0:0.2:20;
H11=diag(ones(1,N_spin-1),1);
H11(1,N_spin)=-1;
H12=H11-transpose(H11);
H12(1,N_spin)=1;
H12(N_spin,1)=-1;
H11=H11+transpose(H11);
H=-J/2*[H11,H12;
    -H12,-H11];
H=H+g*diag([-ones(1,N_spin),ones(1,N_spin)]);

[V,D]=eig(H);

rho0=diag([ones(1,N_spin),-ones(1,N_spin)]);
rho0=rho0/(2*N_spin)^0.5;

rhoV=transpose(V)*rho0*V;

CZX=0;
CZY=0;
L_c1=zeros(size(tlist));
L_c=L_c1;
count=1;
for t=tlist
    rhot=V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*transpose(V);
    rho11=rhot(1:N_spin,1:N_spin);
    rho12=rhot(1:N_spin,N_spin+1:2*N_spin);
    
    CZX=CZX+sum(diag(rho11).^2);
    CZY=CZY+sum(diag(rho11).^2);
    
    M1=(real(rho11+rho12)).^2;
    M2=(real(rho11-rho12)).^2;
    
    for p=-N_spin:N_spin
        L_c(count)=L_c(count)+(abs(p)+1)*sum(abs(diag(rho11,p)).^2);
        L_c(count)=L_c(count)+(abs(p)+1)*sum(abs(diag(rho12,p)).^2);
    end
    L_c1(count)=N_spin*abs(rho11(1,1))^2;
    for p=2:N_spin
        L_c1(count)=L_c1(count)+(N_spin)*(2*p*(abs(rho11(1,p)))^2+2*p*abs(rho12(1,p))^2);
    end
    count=count+1;
end
L_c1=L_c1*2;
figure(1)
plot(tlist,L_c*2)
