%20171018
%Compute M_t of the  integrable model, diagonalizing in
%fermion picture
%H=-gZ-JXX
%Note here H is not the Hamiltonian, it's defined by the propagator in the
%fermion picture: U=e^(iHt). H here is twice the Hamiltonian.
N_spin=10;
J=1;
g=0.02;
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
T=10;
tlist=0:Dt:T;

M_t=tlist;

count=1;
for g=min_g:step_g:max_g
    rhot=V*exp(1i*D*t)*rhoV*diag(exp(diag(-1i*D*t)))*V.';
    M_t=
    count=count+1;
end
L_c=L_c*2;
L_c1=L_c1*2;
figure(1)
plot(x,L_c,x,L_c1)
