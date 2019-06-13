function result = Miecoated_ab2(m1,m2,x,y)

% Basic computation of Mie Coefficients, a_n, b_n, 
% of orders n=1 to nmax, complex refractive index m=m'+im", 
% and size parameters x=k0*a, y=k0*b where k0= wave number 
% in the ambient medium for coated spheres, a=inner radius,
% b=outer radius m1,m2=inner,outer refractive index; 
% p. 183 in Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002

m=m2/m1;
nmax=round(2+y+4*y.^(1/3));
n=(1:nmax); nu = (n+0.5); 
u=m1*x; v=m2*x; w=m2*y;
su=sqrt(0.5*pi*u);sv=sqrt(0.5*pi*v);sw=sqrt(0.5*pi*w);sy=sqrt(0.5*pi*y);
pu=su.*besselj(nu,u); py=sy.*besselj(nu,y);
pv=sv.*besselj(nu,v); pw=sw.*besselj(nu,w);
p1u=[sin(u), pu(1:nmax-1)]; p1y=[sin(y), py(1:nmax-1)];
p1v=[sin(v), pv(1:nmax-1)]; p1w=[sin(w), pw(1:nmax-1)];
ppv=p1v-n.*pv./v; ppw=p1w-n.*pw./w; ppy=p1y-n.*py./y;
chv=-sv.*bessely(nu,v); 
chw=-sw.*bessely(nu,w);
chy=-sy.*bessely(nu,y);
ch1v=[cos(v), chv(1:nmax-1)]; ch1w=[cos(w), chw(1:nmax-1)];
ch1y=[cos(y), chy(1:nmax-1)]; 
gsy=py-i*chy;
gs1y=p1y-i*ch1y; gspy=gs1y-n.*gsy./y;
du=p1u./pu-n./u; dv=p1v./pv-n./v; dw=p1w./pw-n./w;
chpv=ch1v-n.*chv./v; chpw=ch1w-n.*chw./w; 

aan=pv.*(m.*du-dv)./(m.*du.*chv-chpv);
bbn=pv.*(m.*dv-du)./(m.*chpv-du.*chv);

a1=ppw-aan.*chpw;
a2=pw-aan.*chw;
b1=ppw-bbn.*chpw; b2=pw-bbn.*chw; 
an=(py.*a1-m2.*ppy.*a2)./(gsy.*a1-m2.*gspy.*a2);
bn=(m2.*py.*b1-ppy.*b2)./(m2.*gsy.*b1-gspy.*b2);

result=[an; bn];