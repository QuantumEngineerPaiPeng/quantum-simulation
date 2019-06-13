function result = cloudphasefunbeta0(lamda, rc, teta0, factor);

% "beta0" of phase function in Mie Theory (s. Eq. 9 of
% Meador and Weaver, 1980 J. Atm. Sci. 37, pp. 630-643)
% for a cloud with a mod'gamma-size distribution, with
% mode radius rc, wavelength lamda, both in micron, 
% beta0 is a function of incidence angle teta;
% the whole function is computed, and for a given
% value of teta0, the function value is given as beta00
% C. Mätzler, July 2003

dn0=0.7;    % Teta-Versatz der Stützwerte: 0 bis < 1
xc=2*pi*rc/lamda;
nj=ceil((xc+1.5)*2.3);
nsteps=factor*nj;
n2=2*nsteps;
dteta=pi/n2;
nx=(1:n2);
teta=(nx-dn0)*dteta;
u=cos(teta); s=sin(teta);
AA=[];
for j=1:n2,
    ct=u(j); st=s(j);
    cc=cloud_phasefunction(lamda, rc, ct, 5);
    AA=[AA;0.5*dteta*cc(1)*st];
end;
L=[];
for jj=1:nj,
    x0=legendre(jj-1,u);
    x01=x0(1,:);
    L=[L,x01'];    % Legendre Polynom(teta), Order jj-1 in Column jj
    gj(jj)=x01*AA; % gj are the coefficients gl in Mead-Weav
    gj2(jj)=gj(jj)*(2*jj-1);
end;
for j1=1:n2,
    for j2=1:n2,
        ab(j1,j2)=gj2.*L(j1,:)*L(j2,:)';
    end;
end;
ptest=dteta*s*ab(:,1)   % ptest should be 2
nx2=(1:nsteps);
teta2=teta(1:nsteps);
s2=s(1:nsteps);
c2=cos(teta2);
L2=L(1:nsteps,:);
Q=0.5*dteta*s2*L2;
for j=1:nsteps,
    beta0(j)=1-gj2.*L(j,:)*Q';  % beta0(teta0), MeadWeav 1980, Eq. 9
end;
% treatment of boundary values at teta=0 and teta=pi/2:
teta2=[teta2,pi/2]; beta0=[beta0,0.5]; c2=[c2,0]; s2=[s2,1];
nsteps=nsteps+1;
plot(c2,beta0,'k-'), xlabel('mu0'),ylabel('beta0'),
title(sprintf('lamda=%g micron, rc=%g micron',lamda,rc));
% interpolation to find value beta00=beta0(teta0):
t2=min(find(teta2>teta0));
t1=max(find(teta2<=teta0));
if t2==1,
    beta00=beta0(1);
elseif t1==nsteps,
    beta00=0.5;
elseif t1==nsteps-1,
    dt=teta0-pi/2;
    beta01=beta0(nsteps-1);
    beta02=0.5; dteta=pi/2-teta2(nsteps);
    beta00=beta02+dt*(beta02-beta01)/(pi/2-teta2(nsteps-1));
else,
    dt=teta0-teta2(t1);
    beta01=beta0(t1);
    beta02=beta0(t2);
    beta00=beta01+dt*(beta02-beta01)/dteta;
end;
result.beta0=beta00;
result.mu0=c2;
result.beta0fun=beta0;
