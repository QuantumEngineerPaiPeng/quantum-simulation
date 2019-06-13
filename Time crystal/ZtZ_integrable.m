%20180102
%Compute OTOC of the  integrable model, diagonalizing in
%fermion picture
%H=-gZ-JXX
%Note here H is not the Hamiltonian, it's defined by the propagator in the
%fermion picture: U=e^(iHt). H here is twice the Hamiltonian.
N_spin=16;
J=1;
g=0;
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
rho0=rho0/(2*N_spin)^0.5;
rhoV=V.'*rho0*V;

Dt=0.02;
T=10;
tlist=0:Dt:T;

otoc=0;

count=1;
for t=Dt:Dt:T
    rhot=V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*V.';
    otoc=[otoc,-trace((rhot*rho0-rho0*rhot)^2)];
    count=count+1;
end

figure(1)
plot(tlist,otoc)