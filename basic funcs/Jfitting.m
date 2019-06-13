% 20180703
% Get the coupling rate by fitting the dq evolution to a Bessel function
fun=@(a,b,c,gamma,x) b*(besselj(0,a*x)).*exp(-gamma*x)+c;
y=real(B);
x=(1:16)*0.7*96.5;
[f1,f2,f3]=fit(x.',y.',fun,'StartPoint',[0.03,9e3,100,2e-3]);
figure
plot(f1,x,y)
f1