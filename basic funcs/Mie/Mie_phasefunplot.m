function result = mie_phasefunplot(m, x, nsteps)

% Plot of Mie phasefunction for unpolarised radiation 
% with normalisation to 'one' when integrated 
% over all directions/(4*pi), see Chandrasekhar 1960
% Eq. (28), for complex refractive-index ratio m=m'+im", 
% size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% INPUT: 
% m, x: as usual
% nsteps: number of angles to be considered
% C. Mätzler, July 2003.

dteta=pi/(nsteps-1);
m1=real(m); m2=imag(m);
nx=(1:nsteps);
teta=(nx-1).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=Mie_S12(m,x,u);
        SL(j)= real(a(1,j)'*a(1,j));
        SR(j)= real(a(2,j)'*a(2,j));
    end;
tetad=teta*180/pi;
Q=mie(m,x);
Qext=Q(1); Qsca=Q(2); asy=Q(5); w=Qsca/Qext;
p=2*(SL+SR)./Qsca/x.^2;
pdB=10*log10(p);
% Qsca to be exchanged by Qext=Q(1) above if normalisation
% to single-scattering albedo is required.
plot(tetad,pdB,'r-')
title(sprintf('Phase Function: m=%g+%gi, x=%g, w=%g, g=%g',m1,m2,x,w,asy))
xlabel('Scattering Angle (deg)'),ylabel('Phase Function (dB)')
result=[teta;p]'; % Phase function  
