function result = besselplot2(n, m, xmin, dx, nx)

% Computation and plot of Logarithmic Derivatives of Riccati-Bessel Functions 
% of Order n for complex argument z=m*x, 
% input: order n, refractive index m, minimum x value xmin, 
% x interval dx, number of x values nx.
% C. Mätzler, August 2002

m1=real(m);m2=imag(m);
nn=(1:nx);
x=xmin+dx*nn;
nu=n+0.5;
z=m.*x;
sqz= sqrt(0.5*pi*z); 
pz = besselj(nu, z).*sqz;       % Psi_n Function
chz = -bessely(nu, z).*sqz;     % Chi_n Function
p1z= besselj(nu-1, z).*sqz;     % Psi_n-1 Function
ch1z= -bessely(nu-1, z).*sqz;;  % Chi_n-1 Function
dz=p1z./pz-n./z;                % D_n =Psi_n'/Psi_n Function
dcz=ch1z./chz-n./z;             % Dch_n =Chi_n'/Chi_n Function
dz1=real(dz);
dz2=imag(dz);
dcz1=real(dcz);
dcz2=imag(dcz);
r=[dz1;dz2;dcz1;dcz2];

plot(x,r(1:4,:))
legend('real(d_n_p_s(mx))','imag(d_n_p_s(mx))','real(d_n_c_h(mx))','imag(d_n_c_h(mx))')
title(sprintf('Logarithmic Derivatives of Riccati-Bessel Functions of Order n=%g, for m=%g+%gi',n,m1,m2))
xlabel('x')
result=[x;r];