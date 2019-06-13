function result = Mie2_tetascan(eps1, mu1, x, nsteps)

% Computation and plot of Mie Power Scattering function for given 
% complex permittivity and permeability ratios eps1=eps1'+ieps1",
% mu1=mu1'+imu1" between inside and outside of the sphere
% and range of size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002.

nsteps=nsteps;
epsp=real(eps1);epspp=imag(eps1);
mup=real(mu1);mupp=imag(mu1);
nx=(1:nsteps); dteta=pi/(nsteps-1);
teta=(nx-1).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=Mie2_S12(eps1,mu1,x,u);
        SL(j)= real(a(1,j)'*a(1,j))/(pi*x^2);
        SR(j)= real(a(2,j)'*a(2,j))/(pi*x^2);
    end;
y=[teta teta+pi;SL SR(nsteps:-1:1)]'; 

polar(y(:,1),y(:,2))
title(sprintf('Angular pattern for x=%g, eps1=%g+%gi, mu1=%g+%gi',x,epsp,epspp,mup,mupp));
xlabel('Scattering Angle')
result=y; 
