function result = MieRayleighxscan2(m, nsteps, dx, xmax, nmax)

% Computation and plot of Mie Efficiencies for given 
% complex refractive-index ratio m=m'+im" 
% and range of size parameters x=k0*a, 
% starting at x=0 with nsteps increments of dx
% a=sphere radius, using complex Mie coefficients an and bn 
% according to Bohren and Huffman (1983) BEWI:TDD122
% result: m', m", x, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% qratio=qb/qsca and asymmetry parameter (asy=<costeta>).
% nmax is a maximum order of spherical functions used only if
% d is to be computed (line 22)
% C. Mätzler, May 2002.

nx=(1:nsteps)';
x=0.1+(nx-1)*dx;
for j = 1:nsteps,
    a(j,:)=mie(m,x(j));   % Mie Solution
    b(j,:)=mie_1(m,x(j),xmax); % Rayleigh Approximation (1st order)
%   c(j,:)=mie_2(m,x(j),xmax); % Second Order Approximation
%   d(j,:)=mie_nmax(m,x(j),nmax); % nmax Order of Approximation
end;
m1=real(m); m2=imag(m);
%plot(x,a(:,param),'k-',x,b(:,param),'r-',x,c(:,param),'r-.',x,d(:,param),'k:')
%legend('Mie','Rayleigh','2nd Order','n<=nmax')
%title(sprintf('m = %g + %gi,  nmax = %g',m1,m2,nmax))
param=1;
subplot(1,2,1);
loglog(x,a(:,param),'r-',x,b(:,param),'k--')
title(sprintf('n = %g + %gi',m1,m2))
xlabel('x'), ylabel('Qext');
param=4;
subplot(1,2,2);
loglog(x,a(:,param),'r-',x,b(:,param),'k--')
legend('Mie','Rayleigh')
title(sprintf('n = %g + %gi',m1,m2))
xlabel('x'), ylabel('Qb');

%result=a; 