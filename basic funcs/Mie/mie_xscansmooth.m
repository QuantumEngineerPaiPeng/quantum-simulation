function result = mie_xscansmooth(m, nsteps, dx, Dx)

% Computation and plot of Mie Efficiencies for given 
% complex refractive-index ratio m=m'+im" 
% and range of size parameters x=k0*a, 
% starting at x=0 with nsteps increments of dx
% a=sphere radius, using complex Mie coefficients an and bn 
% according to Bohren and Huffman (1983) BEWI:TDD122
% result: m', m", x, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% qratio=qb/qsca and asymmetry parameter (asy=<costeta>).
% C. Mätzler, May 2002.
nsmooth=Dx/dx;
nx=(1:nsteps)';
x=(nx-1)*dx;
for j = 1:nsteps
    a(j,:)=mie(m,x(j));
end;
output_parameters='Qext, Qsca, Qabs, Qb, <costeta>, Qb/Qsca'
for j=1:6,
    b(:,j)=smooth(a(:,j),nsmooth);
end;
b(:,3)=10000*b(:,3);
% plotting the smoothed results
m1=real(m);m2=imag(m);
plot(x,b(:,1:6))
legend('Qext','Qsca','10^4*Qabs','Qb','<costeta>','Qb/Qsca')
title(sprintf('Smoothed (Dx=%g) Mie Efficiencies, m=%g+%gi',Dx,m1,m2))
xlabel('x')

result=b; 