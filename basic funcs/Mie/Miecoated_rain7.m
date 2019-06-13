function result = Miecoated_rain7(fGHz, R, TK, coatmax, nr)

% Extinction, scattering, absorption, backscattering and 
% asymmetric scattering coefficients in 1/km for Marshall-Palmer 
% (MP) drop-size distribution (Sauvageot et al. 1992),
% versus thickness 'coat' of ice-coated water spheres (freezing rain)
% assuming coat=min(coat,radius) for all spheres, using Mie Theory,
% the dielectric model of Liebe et al. 1991 for water and
% of Mätzler (1998) for ice.
% Input:
% coat: thickness of coating in mm, 
% R: rain rate in mm/h, TK: Temp. in K, 
% fmin, fmax: minimum and maximum frequency in GHz
% nr: Number of coat-thickness values 
% C. Mätzler, June 2002.

opt=1;
nsteps=501;                  % number of drop-diameter values optimized for MP 
LA=4.1/R^0.21;               % MP paramter (with size unit in mm)
N0=0.08/10000;               % original MP N0 in 1/mm^4
nx=(1:nsteps)';
c0=299.793;
m1=sqrt(epswater(fGHz, TK)); % refractive index of pure water
m2=sqrt(epsice(fGHz, TK));   % refractive index of pure ice
dD=0.01*R^(1/6)/fGHz^0.05;   % diameter interval optimized for MP
D=(nx-0.5)*dD;               % drop diameter in mm
y=pi*D*fGHz/c0;
dx=2*pi*coatmax*fGHz/c0/(nr-0.99999);
for jr = 1:nr
    x=max(y-dx*(jr-1),0);
    coat=dx*(jr-1)*c0/(2*pi*fGHz);
    sigmag=pi*D.*D/4;        % geometric cross section
    NMP=N0*exp(-LA*D);       % MP distribution
    sn=sigmag.*NMP*1000000;  
    for j = 1:nsteps    
        a(j,:)=Miecoated(m1,m2,x(j),y(j),opt);
    end;
    b(:,1)=D;             b(:,2)=a(:,1).*sn;   
    b(:,3)=a(:,2).*sn;    b(:,4)=a(:,3).*sn;
    b(:,5)=a(:,4).*sn;   b(:,6)=a(:,2).*a(:,5).*sn; 
    gext= sum(b(:,2))*dD;    gsca= sum(b(:,3))*dD;
    gabs= sum(b(:,4))*dD;    gb=   sum(b(:,5))*dD;
    gteta=sum(b(:,6))*dD;
    res(jr,:)=[coat gext gsca gabs gb gteta];
end;
output_parameters='Gext, Gsca, Gabs, Gb, Gsca<costeta>'
semilogy(res(:,1),res(:,2:6))
legend('Gext','Gsca','Gabs','Gb','Gsca<costeta>)')
title(sprintf('Freezing Rain at R=%gmm/h, T=%gK, f=%gGHz',R,TK,fGHz))
xlabel('Ice Coating (mm)');    ylabel('Gi(1/km)')
result=res;
