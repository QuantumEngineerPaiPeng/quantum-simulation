%20180214
%Compute OTOC with one operator being time averaged
N_spin=200;
J=1;
g=3;
%t=0.16*12;

H11=diag(ones(1,N_spin-1),1);
%H11(1,N_spin)=-1;
H12=H11-transpose(H11);
%H12(1,N_spin)=1;
%H12(N_spin,1)=-1;
H11=H11+transpose(H11);
H=-J/2*[H11,H12;
    -H12,-H11];
H=H+g*diag([-ones(1,N_spin),ones(1,N_spin)]);

[V,D]=eig(H);

rho0=diag([ones(1,N_spin),-ones(1,N_spin)]);
%rho0=rho0/(2*N_spin)^0.25;

tlist=12:0.1:50;
rhoV=transpose(V)*rho0*V;

CZX=zeros(size(tlist));
CZY=zeros(size(tlist));

rhobar=rho0;
for t=tlist
    rhobar=rhobar+V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*transpose(V);
end
rhobar=rhobar/(length(tlist)+1);
otoc=-trace((rhobar*rho0-rho0*rhobar)^2)/(2*N_spin)*4;
plot(g,otoc,'k+')