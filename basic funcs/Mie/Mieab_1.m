function result = Mieab_1(m, x)
% Computation of Mie coefficients a_n(z), b_n(z) for order 1,
% complex refractive index m=m'+im", and size parameter x=k0*a<<1, 
% where k0=vacuum wave number, a=sphere radius.
% This is the Rayleigh Approximation
% s. p. 101, 127, 131 in Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, May 2002

z=m*x;
m2=m*m;
k2=(m2-1)/(m2+2);
x2=x*x; x3=x2*x;
a1=-i*2*x3*k2/3;
result=a1;