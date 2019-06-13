function result = Miecoated_ab3(m1,m2,x,y)

% Computation of Mie Coefficients, a_n, b_n, 
% of orders n=1 to nmax, complex refractive index m=m'+im", 
% and size parameters x=k0*a, y=k0*b where k0= wave number 
% in the ambient medium for coated spheres, a=inner radius,
% b=outer radius m1,m2=inner,outer refractive index; 
% p. 183 in Bohren and Huffman (1983) BEWI:TDD122
% but using the bottom equation on p. 483 for chi_prime
% C. Mätzler, August 2002

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
ch1y=[cos(y), chy(1:nmax-1)]; 
gsy=py-i*chy; 
gs1y=p1y-i*ch1y; gspy=gs1y-n.*gsy./y; 
du=p1u./pu-n./u; 
dv=p1v./pv-n./v; 
dw=p1w./pw-n./w;
chpw=chw.*dw-1./pw;
uu=m.*du-dv;
vv=du./m-dv;
pvi=1./pv;
aaa=pv.*uu./(chv.*uu+pvi);
bbb=pv.*vv./(chv.*vv+pvi);
aa1=ppw-aaa.*chpw; aa2=pw-aaa.*chw;
bb1=ppw-bbb.*chpw; bb2=pw-bbb.*chw; 
aa=(py.*aa1-m2.*ppy.*aa2)./(gsy.*aa1-m2.*gspy.*aa2);
bb=(m2.*py.*bb1-ppy.*bb2)./(m2.*gsy.*bb1-gspy.*bb2);

result=[aa; bb];