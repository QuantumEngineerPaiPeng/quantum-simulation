%20171219
%Compute M_t of the  integrable model, diagonalizing in
%fermion picture
%H=-gZ-JXX
%Note here H is not the Hamiltonian, it's defined by the propagator in the
%fermion picture: U=e^(iHt). H here is twice the Hamiltonian.
N_spin=15;
J=1;
g=0.2;
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

Dt=0.1;
T=30;
tlist=0:Dt:T;

M_z=trace(rhoV*rhoV);

count=1;
for t=Dt:Dt:T
    rhot=diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)));
    M_z=[M_z,trace(rhot*rhoV)];
    count=count+1;
end

figure(2)
plot(tlist/2,M_z)