function result = Mie_abcd(m, x)

% Computes a matrix of Mie coefficients, a_n, b_n, c_n, d_n, 
% of orders n=1 to nmax, complex refractive index m=m'+im", 
% and size parameter x=k0*a, where k0= wave number 
% in the ambient medium, a=sphere radius; 
% p. 100, 477 in Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, June 2002

nmax=round(2+x+4*x.^(1/3));
n=(1:nmax); nu = (n+0.5); z=m.*x; m2=m.*m; 
sqx= sqrt(0.5*pi./x); sqz= sqrt(0.5*pi./z);
bx = besselj(nu, x).*sqx;
bz = besselj(nu, z).*sqz;
yx = bessely(nu, x).*sqx;
hx = bx+i*yx;
b1x=[sin(x)/x, bx(1:nmax-1)];
b1z=[sin(z)/z, bz(1:nmax-1)];
y1x=[-cos(x)/x, yx(1:nmax-1)];
h1x= b1x+i*y1x;
ax = x.*b1x-n.*bx;
az = z.*b1z-n.*bz;
ahx= x.*h1x-n.*hx;

an = (m2.*bz.*ax-bx.*az)./(m2.*bz.*ahx-hx.*az);
bn = (bz.*ax-bx.*az)./(bz.*ahx-hx.*az);
cn = (bx.*ahx-hx.*ax)./(bz.*ahx-hx.*az);
dn = m.*(bx.*ahx-hx.*ax)./(m2.*bz.*ahx-hx.*az);
result=[an; bn; cn; dn];