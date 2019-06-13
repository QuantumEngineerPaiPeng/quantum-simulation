function result = Mie2_abs(eps1, mu1, x)

% Computation of the Absorption Efficiency Qabs of a sphere
% for given complex permittivity and permeability ratios,  
% eps1=eps1'+ieps1", mu1=mu1'+imu1", between inside and outside 
% of the sphere for size parameter x=k0*a. 
% Ref. Bohren and Huffman (1983) BEWI:TDD122,
% and my own notes on this topic;
% k0=2*pi./wavelength; % x=k0.*radius;
% Input: eps1, mu1, x
% C. Mätzler, July 2002.

ep2=imag(eps1); mu2=imag(mu1);
nj=100*round(2+x+4*x.^(1/3));
dx=x/nj;
x2=x.*x;
nj1=nj+1;
xj=(0:dx:x);
en=Mie2_Esquare(eps1,mu1,x,nj); % E-field 
en1=0.5*en(nj1).*x2;     % End-Term correction in integral
enx=en*(xj.*xj)'-en1;    % Trapezoidal radial integration
inte=dx.*enx;
Qabse=4.*ep2.*inte./x2
hn=Mie2_Esquare(mu1,eps1,x,nj); % H-field (by duality)
hn1=0.5*hn(nj1).*x2;     % End-Term correction in integral
hnx=hn*(xj.*xj)'-hn1;    % Trapezoidal radial integration
inth=dx.*hnx;
Qabsm=4.*mu2.*inth./x2
Qabs=Qabse+Qabsm;
result=Qabs;
