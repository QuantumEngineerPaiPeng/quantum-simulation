function result = Mie2_xscan(eps1, mu1, nsteps, dx)

% Computation and plot of Mie Efficiencies for given complex 
% permittivity and permeability ratios eps1=eps1'+ieps1",
% mu1=mu1'+imu1" between inside and outside of the sphere
% and range of size parameters x=k0*a, 
% starting at x=0 with nsteps increments of dx
% a=sphere radius, using complex Mie coefficients an and bn 
% according to Bohren and Huffman (1983) BEWI:TDD122
% result: x, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% qratio=qb/qsca and asymmetry parameter (asy=<costeta>).
% C. Mätzler, May 2002.

nx=(1:nsteps)';
x=(nx-1)*dx;
for j = 1:nsteps
    a(j,:)=Mie2(eps1,mu1,x(j));
end;
output_parameters='Real(m), Imag(m), x, Qext, Qsca, Qabs, Qb, <costeta>, Qb/Qsca'

% plotting the results
epsp=real(eps1);epspp=imag(eps1);
mup=real(mu1);mupp=imag(mu1);
plot(x,a(:,1:6))
legend('Qext','Qsca','Qabs','Qb','<costeta>','Qb/Qsca')
title(sprintf('Mie Efficiencies for eps1=%g+%gi, mu1=%g+%gi',epsp,epspp,mup,mupp))
xlabel('x')

result=[x a]; 