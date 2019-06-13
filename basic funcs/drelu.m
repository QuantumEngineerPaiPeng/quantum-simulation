function y=drelu(dx,x)
y=zeros(size(x));
y(x>0)=dx(x>0);