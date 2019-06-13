%20171016
%Compute correlation length of transverse Ising model, with density matrix
%analytically obtained
N_spin=20;
U=zeros(N_spin,N_spin);
V=U;
F=U;
J=1;
g=0.1;
%t=0.16*12;


plist=-(N_spin-1)/2:1:(N_spin-1)/2;
plist=plist*2*pi/N_spin;

count=1;
t=1;
p=1;
r=1;
ur=0;
vr=0;
rho0r=0;
rho2r=0;
for p=1:N_spin
    w_p=(J^2+g^2+2*J*g*cos(plist(p)))^0.5;
    sin_theta_p=J*sin(plist(p))/w_p;
    cos_theta_p=(J*cos(plist(p))+g)/w_p;
    ur=ur+1/N_spin*exp(-1i*r*plist(p))*(cos(w_p*t)-1i*cos_theta_p*sin(w_p*t));
    vr=vr+1/N_spin*exp(1i*r*plist(p))*(-sin_theta_p*sin(w_p*t));
    rho0r=rho0r+1/N_spin*exp(1i*r*plist(p))*(cos(w_p*t)^2+(cos_theta_p^2-sin_theta_p^2)*sin(w_p*t)^2);
    rho2r=rho2r+2/N_spin*exp(-1i*r*plist(p))*(cos(w_p*t)+1i*cos_theta_p*sin(w_p*t))*sin_theta_p*sin(w_p*t);
end
rho2r