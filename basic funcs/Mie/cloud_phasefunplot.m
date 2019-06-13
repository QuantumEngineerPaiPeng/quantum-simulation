function result = cloud_phasefunplot(lamda, rc, nsteps)

% Plot of Phase Function pm (unpolarised)
% of clouds with normalisation to 'one' when integrated 
% over all directions/(4*pi), see Chandrasekhar 1960
% Eq. (28), 
% Input: wavelength lamda, mode radius rc, both in micron, 
% nsteps=number of scattering angles  
% Output: nstep pairs of: scattering angle, phase function  
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2003

dteta=pi/(nsteps-1);
teta=(0:dteta:pi)';
p=[];
for j=1:nsteps,
    mu=cos(teta(j));
    c=cloud_phasefunction(lamda, rc, mu, 17);
    % last number is number of drop radii in size distribution
    p=[p; c(1)];
end;
semilogy(teta*180/pi,p,'r-')
title(sprintf('Phase Function of Cloud with rc=%g micron, at lamda=%g micron ',rc,lamda))
xlabel('Scattering Angle (deg)'),ylabel('Phase Function')
% result=[teta, p];