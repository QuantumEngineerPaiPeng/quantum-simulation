function result = mie_beamefficiency(m, x, tetalimit, nsteps)

% Cumulative power (corresponds to the beam efficiency in 
% antenna theory) scattered within a maximum scattering angle 
% (variable from 0 to pi) for Mie Scattering with
% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% INPUT: 
% m, x: as usual
% tetalimit: maximum angle to be considered, in radians
% nsteps: number of angles to be considered
% C. Mätzler, June 2003.

dteta=tetalimit/nsteps;
xmin=dteta/pi;
m1=real(m); m2=imag(m);
nx=(1:nsteps); dteta=tetalimit/nsteps;
teta=(nx-0.5).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=Mie_S12(m,x,u);
        SL(j)= real(a(1,j)'*a(1,j))/(pi*x^2);
        SR(j)= real(a(2,j)'*a(2,j))/(pi*x^2);
    end;
st=2*pi*sin(teta);
SSL=st.*SL;
SSR=st.*SR;
tetad=teta*180/pi;
Q=mie(m,x);
Qsca=Q(2);
z=0.5*dteta*cumsum(SSL+SSR)/Qsca;
semilogx(tetad/180,z,'r-')
title(sprintf('Cumulative Fraction of Scattered Power: m=%g+%gi, x=%g',m1,m2,x))
xlabel('Maximum Scattering Angle/180°'),
axis([xmin, tetalimit/pi, 0, 1.1]);
result=[teta; SSL; SSR; z]';
