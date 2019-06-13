function result = besselplot3(n, m, xmin, dx, nx)

% Computation and plot of Inverse Products of Riccati-Bessel Functions 
% of Order n for complex argument z=m*x, used in Mie Theory, 
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
pc=1./(pz.*chz);
dz1=real(pc);
dz2=imag(pc);
r=[dz1;dz2];

plot(x,r(1:2,:))
legend('real(1/(Psi_n(mx)*Chi_n(mx)))','imag(1/(Psi_n(mx)*Chi_n(mx)))')
title(sprintf('Inverse product of Riccati-Bessel Functions of Order n=%g, for m=%g+%gi',n,m1,m2))
xlabel('x')
result=[x;r];