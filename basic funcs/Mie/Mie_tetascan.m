function result = Mie_tetascan(m, x, nsteps, type)

% Computation and plot of Mie Power Scattering function for given 
% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% type ='pol' for polar diagram, else cartesian
% C. Mätzler, May 2002.

nsteps=nsteps;
m1=real(m); m2=imag(m);
nx=(1:nsteps); dteta=pi/(nsteps-1);
teta=(nx-1).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=Mie_S12(m,x,u);
        SL(j)= real(a(1,j)'*a(1,j))/(pi*x^2);
        SR(j)= real(a(2,j)'*a(2,j))/(pi*x^2);
    end;
y=[teta teta+pi;SL SR(nsteps:-1:1)]'; 
tetad=teta*180/pi;
if type=='pol'
    polar(y(:,1),y(:,2))
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x));
    xlabel('Scattering Angle')
else
    semilogy(tetad,SR,'r-',tetad,SL,'b--')
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x))
    xlabel('Scattering Angle')
    legend('SR:','SL:')
end;
result=y; 
