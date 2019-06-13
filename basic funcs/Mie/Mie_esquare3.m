function result = Mie_Esquare3(m, x, n)

% Test of term-wise equality of absorption efficiency
% as computed by two different methods
% n is the term number of the vector-wave mode,
% x size parameter, m refractive index of sphere, 
% Qabsa=abs.efficiency due to an-term
% Qabsb=abs.efficiency due to bn-term
% Qabsc=abs.efficiency due to cn-term
% Qabsd=abs.efficiency due to dn-term
% Sept. 13, 2002

nj=100*round(2+x+4*x.^(1/3))+100;
nu =(n+0.5); 
m1=real(m); m2=imag(m);e2=imag(m.*m);
abcd=Mie_ab(m,x);
an=abcd(1,n);bn=abcd(2,n);
an2=abs(an).^2;
bn2=abs(bn).^2;

abcd=Mie_cd(m,x);
cn=abcd(1,n);dn=abcd(2,n);
cn2=abs(cn).^2;
dn2=abs(dn).^2;
dx=x/nj;
cs=[];ds=[];
for j=1:nj,
    xj=dx.*j;
    z=m.*xj;
    sqz= sqrt(0.5*pi./z);
    bz = besselj(nu, z).*sqz;      % This is j_n(z)
    bz2=(abs(bz)).^2;
    b1z = besselj(nu-1, z).*sqz;   % This is j_n-1(z)
    az = b1z-n.*bz./z;
    az2=(abs(az)).^2;
    z2=(abs(z)).^2;
    n1 =n*(n+1);
    n2 =2*(2*n+1);
    mn=real(bz2.*n2);
    nn1=az2;
    nn2=bz2.*n1./z2;
    nn=n2.*real(nn1+nn2);
    cs=[cs 0.25*(cn2*mn)]; 
    ds=[ds 0.25*(dn2*nn)]; 
end;
xxj=[0:dx:x];
ccn=[cs(1) cs];
ddn=[ds(1) ds];
ef=[ccn;ddn];
plot(xxj,ef(1:2,:));
title(sprintf('Sq. Ampl. Field in Sphere of cn and dn Modes, n=%g, m=%g+%gi, x=%g',n,m1,m2,x))
legend('cn Term','dn Term')
xlabel('r k')

x2=x.*x;
nj1=nj+1;
cn1=0.5*ccn(nj1).*x2;     % End-Term correction in integral
cnx=ccn*(xxj.*xxj)'-cn1;    % Trapezoidal radial integration
intc=dx.*cnx;
Qabsc=4.*e2.*intc./x2;
dn1=0.5*ddn(nj1).*x2;     % End-Term correction in integral
dnx=ddn*(xxj.*xxj)'-dn1;    % Trapezoidal radial integration
intd=dx.*dnx;
Qabsd=4.*e2.*intd./x2;

Qabsb=(n2/x2)*(real(bn)-bn2);
Qabsa=(n2/x2)*(real(an)-an2);
result=[Qabsa;Qabsb;Qabsc;Qabsd];
