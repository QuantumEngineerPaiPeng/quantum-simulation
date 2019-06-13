%20171018
%compute the analytically obtained correlation length for transverse Ising
%model when g<<J.
function y=left_limit(J,g,t)
a=g/J;
x=2*J*a*t;
J0=besselj(0,x);
J1=besselj(1,x);
y=1/6*(3+(3+4*x.^2).*J0.^2-4*x.*J0.*J1+(1+4*x.^2).*J1.^2)...
    -cos(2*J*t).*(J0+besselj(2,x))+1;
end