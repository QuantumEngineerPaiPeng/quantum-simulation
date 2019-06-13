function result = besselplot5(n, m, xmin, dx, nx)

% Computation and plot of absolute values of 
% Riccati-Bessel Functions of Order n
% for complex argument z=m*x, used in Mie Theory. 
% input: order n, refractive index m, minimum x value xmin, 
% x interval dx, number of x values nx.
% C. Mätzler, August 2002

m1=real(m); m2=imag(m);
nn=(1:nx)
x=xmin+dx*nn;
nu=n+0.5;
z=m.*x;
sqx= sqrt(0.5*pi*z); 
psx = besselj(nu, z).*sqx
chx = -bessely(nu, z).*sqx;
a1=abs(psx);
a2=abs(chx);
r=[a1;1./a2];
semilogy(x,r(1:2,:))
    legend('abs(psi_n(mx))','1/abs(chi_n(mx))')
    title(sprintf('Riccati-Bessel Functions of Order n=%g, for m=%g+%gi',n,m1,m2))
    xlabel('x')
end;
result=[a1(1),a1(nx);a2(1),a2(nx)];