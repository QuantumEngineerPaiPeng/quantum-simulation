%20160719
%Distribution of x of a harmonic oscillator under Boltzmann distribution
min=1e-10;
max=5;
N=100;
step=(max-min)/N;

x=min:step:max;

y=zeros(N+1,1);

trun=100;

for p=1:N
    y(p)=4*dblquad('ReA',x(p),trun,0,trun);
end

plot(x,y)