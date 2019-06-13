function result = Mie_abs(m, x)

% Computation of the Absorption Efficiency Qabs
% of a sphere of size parameter x, 
% complex refractive index m=m'+im", 
% based on nj internal radial electric field values
% to be computed with Mie_Esquare(nj,m,x)
% Ref. Bohren and Huffman (1983) BEWI:TDD122,
% and my own notes on this topic;
% k0=2*pi./wavelength;
% x=k0.*radius;
% C. Mätzler, May 2002, revised July 2002.

nj=100*round(2+x+4*x.^(1/3))+300;
e2=imag(m.*m);
dx=x/nj;
x2=x.*x;
nj1=nj+1;
xj=(0:dx:x);
en=Mie_Esquare(m,x,nj);
en1=0.5*en(nj1).*x2;     % End-Term correction in integral
enx=en*(xj.*xj)'-en1;    % Trapezoidal radial integration
inte=dx.*enx;
Qabs=4.*e2.*inte./x2;
result=Qabs;
