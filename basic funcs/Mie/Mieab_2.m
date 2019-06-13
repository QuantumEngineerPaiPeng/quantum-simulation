function result = Mieab_2(m, x)

% Low-frequency approximation of Mie Coeffs. a_n(z), b_n(z), order n=1,2
% complex refractive index m=m'+im", and size parameter x=k0*a<<1, 
% where k0=vacuum wave number, a=sphere radius,
% p. 101, 127, 131 in Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, May 2002

z=m*x;
m2=m*m;
k2=(m2-1)/(m2+2);
x2=x*x; x3=x2*x;
a_n(1)=-i*2*x3*k2*(1/3+x2*(m2-2)/(m2+2)/5)+(2*x3*k2/3)^2;
b_n(1)=-i*x2*x3*(m2-1)/45;
a_n(2)=-i*x2*x3*(m2-1)/(2*m2+3)/15;
b_n(2)=0;
result=[a_n; b_n];