function result = Miecoated_rain6(fGHz, R, TK)

% Extinction, scattering, absorption, backscattering and 
% asymmetric scattering coefficients in 1/km for Marshall-Palmer 
% (MP) drop-size distribution (Sauvageot et al. 1992),
% versus thickness 'coat' of water-coated ice spheres (melting ice)
% assuming coat=min(coat,radius) for all spheres, using Mie Theory,
% the dielectric model of Liebe et al. 1991 for water and
% of Mätzler (1998) for ice.
% Input:
% fGHz: frequency in GHz, R: rain rate in mm/h, TK: Temp. in K
% C. Mätzler, July 2002.

opt=1;
nsteps=201;                  % number of drop-diameter values optimized for MP 
LA=4.1/R^0.21;               % MP paramter (with size unit in mm)
N0=0.08/10000;               % original MP N0 in 1/mm^4
nx=(1:nsteps)';
c0=299.793;
m2=sqrt(epswater(fGHz, TK)); % refractive index of pure water
m1=sqrt(epsice(fGHz, TK));   % refractive index of pure ice
dD=0.025*R^(1/6)/fGHz^0.05;   % diameter interval optimized for MP
D=(nx-0.5)*dD;               % drop diameter in mm
y=pi*D*fGHz/c0;
coa=[0.,0.000001,0.000003,0.00001,0.00003,0.0001,0.0003,0.001,0.002,0.004,0.008,0.012,0.02,0.03,0.05,0.08,0.12,0.18,0.27,0.38,0.60,0.75,1];
for jr = 1:23
    dx=2*pi*coa(jr)*fGHz/c0;
    x=max(y-dx,0);
    coat=dx*c0/(2*pi*fGHz);
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
loglog(res(:,1),res(:,2:6))
legend('Gext','Gsca','Gabs','Gb','Gsca<costeta>)')
title(sprintf('Melting ice rain at R=%gmm/h, T=%gK, f=%gGHz',R,TK,fGHz))
xlabel('Water Coating (mm)');    ylabel('Gi(1/km)')
result=res;
