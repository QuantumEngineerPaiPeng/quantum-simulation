function result = Mie_Esquare2(m, x, nj)

% Computation of nj+1 equally spaced values within (0,x)
% of the mean-absolute-square internal 
% electric field of a sphere of size parameter x, 
% complex refractive index m=m'+im", 
% where the averaging is done over teta and phi,
% with unit-amplitude incident field;
% Ref. Bohren and Huffman (1983) BEWI:TDD122,
% and my own notes on this topic;
% k0=2*pi./wavelength;
% x=k0.*radius;
% C. Mätzler, May 2002, revised July 2002

nmax=round(2+x+4*x^(1/3));
n=(1:nmax);  nu =(n+0.5); 
m1=real(m); m2=imag(m);
abcd=Mie_cd(m,x);
cn=abcd(1,:);dn=abcd(2,:);
cn2=abs(cn).^2;
dn2=abs(dn).^2;
dx=x/nj;
es=[];
for j=1:nj,
    xj=dx.*j;
    z=m.*xj;
    sqz= sqrt(0.5*pi./z);
    bz = besselj(nu, z).*sqz;      % This is jn(z)
    bz2=(abs(bz)).^2;
    b1z=[sin(z)/z, bz(1:nmax-1)];  % Note that sin(z)/z=j0(z)
    az = b1z-n.*bz./z;
    az2=(abs(az)).^2;
    z2=(abs(z)).^2;
    n1 =n.*(n+1);
    n2 =2.*(2.*n+1);
    mn=real(bz2.*n2);
    nn1=az2;
    nn2=bz2.*n1./z2;
    nn=n2.*real(nn1+nn2);
    es=[es 0.25*(cn2*mn'+dn2*nn')];
end;
xxj=[0:dx:xj]; een=[es(1) es];

plot(xxj,een);

legend('Radial Dependence of (abs(E))^2')
title(sprintf('Squared Amplitude E Field in a Sphere, m=%g+%gi, x=%g',m1,m2,x))
xlabel('r k')

result=een;    
