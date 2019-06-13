function result = Miecoated_prscan(m1,m2,y,nsteps)

% Computation and plot of Mie Efficiencies of coated
% spheres for given complex refractive-index ratios
% m1,2=m1',2'+im1",2" in kernel and coating, resp.
% and size parameters x=k0*a, y=k0*b, vs. x/y=1-wr (see below)
% with nsteps increments 
% a,b=inner,outer sphere radius, using complex Mie coefficients 
% an and bn, according to Bohren and Huffman (1983) BEWI:TDD122
% result: m1,m2,x,y,wr, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% qratio=qb/qsca, asymmetry parameter (asy=<costeta>),
% volume fraction of coating w= 1-(a/b)^3 = 1-(x/y)^3
% C. Mätzler, July 2002.

m1p=real(m1); m1pp=imag(m1);
m2p=real(m2); m2pp=imag(m2);
nx=(1:nsteps)';
pr=(nx-0.5)/nsteps;               % note that pr=a/b=(1-w)^(1/3)
x=y*pr;
for j = 1:nsteps
    a(j,:)=Miecoated(m1,m2,x(j),y,1);
end;
output_parameters='m1(real,imag), m2(real, imag), x, y, Qext, Qsca, Qabs, Qb, <costeta>, Qb/Qsca'

% plotting the results
plot(pr,a(:,1:5))
legend('Qext','Qsca','Qabs','Qb','<costeta>')
title(sprintf('Mie Efficiencies of coated sphere, y=%g, m1=%g+%gi, m2=%g+%gi ',y,m1p,m1pp,m2p,m2pp))
xlabel('a/b')

result=a; 