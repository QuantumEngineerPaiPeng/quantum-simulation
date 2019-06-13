function result = Mie2_S12(eps1, mu1, x, u)

% Computation of Mie Scattering functions S1 and S2 for
% complex permittivity and permeability ratios eps1=eps1'+ieps1",
% mu1=mu1'+imu1" between inside and outside of the sphere
% and range of size parameters x=k0*a, a=sphere radius
% according to Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002.

nmax=round(2+x+4*x^(1/3));
ab=Mie2_ab(eps1,mu1,x);
an=ab(1,:);
bn=ab(2,:);

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