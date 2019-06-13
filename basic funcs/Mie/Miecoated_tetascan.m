function result = Miecoated_tetascan(m1, m2, x, y, nsteps)

% Computation and plot of Mie Power Scattering function for given 
% complex refractive-index ratios m1,2=m1,2'+im1,2", size parameters 
% x=k0*a, y=k0*b, of coated sphere
% according to Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2002.

nsteps=nsteps;
m1p=real(m1); m1pp=imag(m1);
m2p=real(m2); m2pp=imag(m2);
nx=(1:nsteps); dteta=pi/(nsteps-1);
teta=(nx-1).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=Miecoated_S12(m1,m2,x,y,u);
        SL(j)= real(a(1,j)'*a(1,j))/(pi*y^2);
        SR(j)= real(a(2,j)'*a(2,j))/(pi*y^2);
    end;
scan=[teta teta+pi;SL SR(nsteps:-1:1)]'; 

polar(scan(:,1),scan(:,2)),
title(sprintf('Angular scattering of coated sphere: m1=%g+%gi, m2=%g+%gi, x=%g, y=%g',m1p,m1pp,m2p,m2pp,x,y)),
xlabel('Scattering Angle');
result=y; 
