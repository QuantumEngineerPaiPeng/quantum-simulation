function result = Mie_cd(m, x)

% Computes a matrix of Mie coefficients, c_n, d_n, 
% of orders n=1 to nmax, complex refractive index m=m'+im", 
% and size parameter x=k0*a, where k0= wave number 
% in the ambient medium, a=sphere radius; 
% p. 100, 477 in Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002

z=m.*x;
nmax=round(2+x+4*x.^(1/3));
nmx=round(max(nmax,abs(z))+16);
n=(1:nmax); nu = (n+0.5); msq=m.*m; 
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
nenn1=msq.*ahx-hx.*cnn;
nenn2=ahx-hx.*cnn;
cn = bz.*(bx.*ahx-hx.*ax)./nenn2;
dn = bz.*m.*(bx.*ahx-hx.*ax)./nenn1;
result=[cn; dn];