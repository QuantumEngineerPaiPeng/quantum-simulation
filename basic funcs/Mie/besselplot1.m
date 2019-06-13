function result = besselplot1(n, m, xmin, dx, nx)

% Computation and plot of Riccati-Bessel Functions of Order n
% for complex argument z=m*x, used in Mie Theory. 
% input: order n, refractive index m, minimum x value xmin, 
% x interval dx, number of x values nx.
% C. Mätzler, August 2002

m1=real(m); m2=imag(m);
nn=(1:nx);
x=xmin+dx*nn;
nu=n+0.5;
z=m.*x;
sqx= sqrt(0.5*pi*z); 
psx = besselj(nu, z).*sqx;
chx = -bessely(nu, z).*sqx;
dpic=psx-i*chx;
r=[real(psx);-imag(psx);real(chx);imag(chx);real(dpic);imag(dpic)];
if abs(m2)<0.1
    plot(x,r(1:4,:))
    legend('real(psi_n(mx))','-imag(psi_n(mx))','real(chi_n(mx))','imag(chi_n(mx))')
    title(sprintf('Riccati-Bessel Functions of Order n=%g, for m=%g+%gi',n,m1,m2))
    xlabel('x')
else
    semilogy(x,r(1:6,:))
    legend('real(psi_n(mx))','-imag(psi_n(mx))','real(chi_n(mx))','imag(chi_n(mx))','real(dpic(mx))','imag(dpic(mx))')
    title(sprintf('Riccati-Bessel Functions of Order n=%g, for m=%g+%gi',n,m1,m2))
    xlabel('x')
end;
result=[x;r];