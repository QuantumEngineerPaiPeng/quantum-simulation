function result = Miecoated_S12(m1, m2, x, y, u)

% Coated-sphere Mie-Scattering functions S1 and S2
% for complex refractive index ratios m1,2=m1,2'+im1,2", 
% size parameters x=k0*a, y=k0*b, and u=cos(scattering angle),
% where k0=vacuum wave number, a=sphere radius;
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002

nmax=round(2+y+4*y^(1/3));
abcd=Miecoated_ab1(m1,m2,x,y);
an=abcd(1,:);
bn=abcd(2,:);

pt=Mie_pt(u,nmax);
pin =pt(1,:);
tin =pt(2,:);

n=(1:nmax);
n2=(2*n+1)./(n.*(n+1));
pin=n2.*pin;
tin=n2.*tin;
S1=(an*pin'+bn*tin');
S2=(an*tin'+bn*pin');
    
result=[S1;S2];