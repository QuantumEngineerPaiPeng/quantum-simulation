function result = Miecoated_rain5(coat, R, TK, fmin, fmax, nfreq)

% Freezing rain
% Extinction, scattering, absorption, backscattering and 
% asymmetric scattering coefficients in 1/km for Marshall-Palmer 
% (MP) drop-size distribution (Sauvageot et al. 1992), versus 
% frequency of ice-coated water spheres (freezing rain)
% with constant thickness 'coat' of the coating, using Mie Theory,
% the dielectric model of Liebe et al. 1991 for water and
% of Mätzler (1998) for ice.
% Input:
% coat: thickness of coating in mm, 
% R: rain rate in mm/h, TK: Temp. in K, 
% fmin, fmax: minimum and maximum frequency in GHz
% nfreq: Number of frequencies 
% C. Mätzler, June 2002.

opt=1
nsteps=501;                  % number of drop-diameter values
N0=0.08/10000;               % original MP N0 in 1/mm^4
fact=(fmax/fmin)^(1/(nfreq-0.99999));
fGHz=fmin/fact;
nx=(1:nsteps)';
c0=299.793;
for jr = 1:nfreq
    fGHz=fGHz*fact;
    m1=sqrt(epswater(fGHz, TK));
    m2=sqrt(epsice(fGHz, TK));
    dD=0.01*R^(1/6)/fGHz^0.05;
    D=(nx-0.5)*dD;
    y=pi*D*fGHz/c0;
    dxcoat=2*pi*coat.*fGHz/c0;
    x=max(0,y-dxcoat);
    sigmag=pi*D.*D/4;
    LA=4.1/R^0.21;
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
    res(jr,:)=[fGHz gext gsca gabs gb gteta];
end;
output_parameters='Gext, Gsca, Gabs, Gb, Gsca<costeta>'
plot(res(:,1),res(:,2:6))
legend('Gext','Gsca','Gabs','Gb','Gsca<costeta>')
title(sprintf('Ice-Coated Rain Coefficients vs. Frequency at R=%gmm/h, T=%gK, coat=%gmm',R,TK,coat))
xlabel('f (GHz)');    ylabel('Gi(1/km)')
result=res;
