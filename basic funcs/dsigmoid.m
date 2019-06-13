function y=dsigmoid(dx,x)
s=sigmoid(x);
y=dx.*s.*(1-s);