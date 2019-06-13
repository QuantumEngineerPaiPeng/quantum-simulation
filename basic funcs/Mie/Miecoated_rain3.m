function result = Miecoated_rain3(fGHz, TK, nrain, pam, coat)

% Melting rain
% Extinction, scattering, absorption, backscattering and 
% asymmetric scattering coefficients in 1/km versus rain rate, 
% for Marshall-Palmer (MP) drop-size distribution 
% see Sauvageot et al. (1992),
% using Mie Theory, and Liebe '91 dielectric model. Input:
% fGHz: frequency in GHz, TK: Temp. in K, 
% nrain: Number of rain rates between Rmin=0.1 and Rmax=100mm/h
% pam: 0 if no costeta data to be given, 1 if they are needed
% coat: thickness of water coating of ice sphere in mm
% C. Mätzler, July 2002.

opt=1;     
Rmin=0.1;   
nsteps=201;
m1=sqrt(epsice(fGHz, TK));
m2=sqrt(epswater(fGHz, TK));
N0=0.08/10000;              % original MP N0 in 1/mm^4
fact=1000^(1/(nrain-0.99999));
R=Rmin/fact;
nx=(1:nsteps)';
c0=299.793;
for jr = 1:nrain
    R=R*fact;
    dD=0.025*R^(1/6)/fGHz^0.05;
    D=(nx-1)*dD;
    y=pi*D*fGHz/c0;
    x=max(0,pi*(D-coat)*fGHz/c0);
    sigmag=pi*D.*D/4;
    LA=4.1/R^0.21;
    NMP=N0*exp(-LA*D);
    sn=sigmag.*NMP*1000000;
    for j = 1:nsteps    
        a(j,:)=Miecoated(m1,m2,x(j),y(j),opt);
    end;
    b(:,1)=D;
    b(:,2)=a(:,1).*sn;
    b(:,3)=a(:,2).*sn;
    b(:,4)=a(:,3).*sn;
    b(:,5)=a(:,4).*sn;
    b(:,6)=a(:,2).*a(:,5).*sn; 
    gext= sum(b(:,2))*dD;
    gsca= sum(b(:,3))*dD;
    gabs= sum(b(:,4))*dD;
    gb=   sum(b(:,5))*dD;
    gteta=sum(b(:,6))*dD;
    res(jr,:)=[R gext gsca gabs gb gteta];
end;
    if pam==0
        output_parameters='Gext, Gsca, Gabs, Gb';
        loglog(res(:,1),res(:,2:5))
        legend('Gext','Gsca','Gabs','Gb')
        title(sprintf('Propagation Coefficients Versus Melting-Rain Rate at f=%gGHz, T=%gK, coat=%gmm',fGHz,TK,coat))
        xlabel('R (mm/h)');    ylabel('Gi(1/km)')
    elseif pam==1
        output_parameters='Gext, Gsca, Gabs, Gb, Gsca*<costeta>';
        loglog(res(:,1),res(:,2:6))
        legend('Gext','Gsca','Gabs','Gb','Gsca*<costeta>')
        title(sprintf('Propagation Coefficients Versus Melting-Rain Rate at f=%gGHz, T=%gK, coat=%gmm',fGHz,TK,coat))
        xlabel('R (mm/h)');    ylabel('Gi(1/km)')
    end;
result=res;
