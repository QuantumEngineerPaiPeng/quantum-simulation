function result = besselriccati(n, m, x)

% Computation of Riccati-Bessel Functions of Order n
% for complex argument z=m*x, used in Mie Theory. 
% input: order n, refractive index m, x value, 
% C. Mätzler, August 2002

m1=real(m); m2=imag(m);
nu=n+0.5;
z=m.*x;
sqx= sqrt(0.5*pi*z); 
psz = besselj(nu, z).*sqx;
chz = -bessely(nu, z).*sqx;
dpic=psz-i*chz;
z2=0.5*z.*z;
n2=2*n;
A1=4./((n2-1).*(n2+3));
lez=(2*n+1)./z./(1+z2.*A1); % Low-frequency approximation of Ez
ez=1./(chz.*psz);
fz=psz./chz/i-1;
result=[psz;chz;dpic;ez;fz;];