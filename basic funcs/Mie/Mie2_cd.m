function result = Mie2_cd(eps1, mu1, x)

% Computes a matrix of Mie coefficients, c_n, d_n, of orders
% n=1 to nmax, for given complex permittivity and permeability
% ratios eps1=eps1'+ieps1", mu1=mu1'+imu1" between inside and outside 
% of the sphere and size parameter x=k0*a 
% where k0= wave number in ambient medium for spheres of radius a;
% p. 100, 477 in Bohren and Huffman (1983), BEWI:TDD122
% C. Mätzler, July 2002

msq=eps1.*mu1;
m=sqrt(msq);            % refractive index ratio
z1=sqrt(mu1./eps1);           % impedance ratio
nmax=round(2+x+4*x.^(1/3));
nmx=round(max(nmax,abs(m))+16);
n=(1:nmax); nu = (n+0.5); z=m.*x;
cnx(nmx)=0+0i;
for j=nmx:-1:2
    cnx(j-1)=j-z.*z./(cnx(j)+j);
end;
cnn=cnx(n);
sqx= sqrt(0.5*pi./x); sqz= sqrt(0.5*pi./z);
bx = besselj(nu, x).*sqx;
bz = 1./besselj(nu, z)./sqz;
yx = bessely(nu, x).*sqx;
hx = bx+i*yx;
b1x=[sin(x)/x, bx(1:nmax-1)];
y1x=[-cos(x)/x, yx(1:nmax-1)];
h1x= b1x+i*y1x;
ax = x.*b1x-n.*bx;
ahx= x.*h1x-n.*hx;
nenn1=msq.*ahx-mu1.*hx.*cnn;
nenn2=mu1.*ahx-hx.*cnn;
cn = mu1.*bz.*(bx.*ahx-hx.*ax)./nenn2;
dn = mu1.*m.*bz.*(bx.*ahx-hx.*ax)./nenn1;
result=[cn; dn];