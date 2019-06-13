function result = besselplot4(n, m, xmin, dx, nx)

% Computation and plot of difference of Logarithmic Derivatives 
% of Riccati-Bessel Functions Psi_n of Order n for complex 
% argument z=m*x, computed in 2 different ways: dz is computed
% from the built-in Matlab Bessel Functions, and dn from the 
% recurrence relation (4.89) of BH.
% input: order n, refractive index m, minimum x value xmin, 
% x interval dx, number of x values nx.
% C. Mätzler, August 2002

xmax=xmin+dx*nx;
m1=real(m);m2=imag(m);
nn=(1:nx);
x=xmin+dx*nn;
nu=n+0.5;
z=m.*x;
zmax=m.*xmax;
sqz= sqrt(0.5*pi*z); 
pz = besselj(nu, z).*sqz;       % Psi_n Function
p1z= besselj(nu-1, z).*sqz;     % Psi_n-1 Function
dz=p1z./pz-n./z;                % D_n =Psi_n'/Psi_n Function

nmx= round(max(2+xmax+4*xmax.^(1/3),abs(zmax))+16);
dnx(1:nmx,nn)=0;
dnx(nmx,nn)=0+0i;
for j=nmx:-1:n+1      % Computation of Dn(z) according to (4.89) of B+H (1983)
    dnx(j-1,nn)=j./z-1./(dnx(j,nn)+j./z);
end;
dn=dnx(n,nn);

dz1=real(dz-dn);
dz2=imag(dz-dn);

r=[dz1;dz2];

plot(x,r(1:2,:))
legend('real(deltad_n_p_s(mx))','imag(deltad_n_p_s(mx))')
title(sprintf('Differences of Log. Derivatives n=%g, for m=%g+%gi',n,m1,m2))
xlabel('x')
% result=[dz1;dz2];
result=[x;dn;log10(pz)];

