function result = phasefunav_phi(m,x,factor);

% Legendre coefficients gi of phase function p(costeta) 
% in Mie Theory (e.g. gl in
% Meador and Weaver, 1980 J. Atm. Sci. 37, pp. 630-643)

% Input: 
% m: refractive index, x: size parameter
% nj: number of gi values (use about nj=2.5*(x+1))
% nsteps: number of angular values between 0 and pi
% (take about nsteps = 20*nj)
% C. Mätzler, July 2003

nj=ceil((x+1.5)*2.3);
nsteps=factor*nj
dteta=pi/nsteps;
m1=real(m); m2=imag(m);
nx=(1:nsteps);
teta=(nx-0.75)*dteta; u=cos(teta); s=sin(teta);
    for j = 1:nsteps, 
        a(:,j)=Mie_S12(m,x,u(j));
        SL(j)= real(a(1,j)'*a(1,j));
        SR(j)= real(a(2,j)'*a(2,j));
        A(j)=(SL(j)+SR(j)).*s(j);
    end;
Q=mie(m,x);
Qext=Q(1); Qsca=Q(2); asy=Q(5); w=Qsca/Qext;
AA=dteta*A./Qsca/x.^2;

L=[];
for jj=1:nj,
    x0=legendre(jj-1,u);
    x=x0(1,:);
    L=[L,x'];    % Legendre Polynom(teta), Order jj-1 in Column jj
    gj(jj)=x*AA'; % gj are the coefficients gl in Mead-Weav
    gj2(jj)=gj(jj)*(2*jj-1);
end;
for j1=1:nsteps,
    for j2=1:nsteps,
        ab(j1,j2)=gj2.*L(j1,:)*L(j2,:)';
    end;
end;
%test: 
ptest=dteta*s*ab(:,1)
% ab is the matrix of (symmetrical) p(mu0,mus) values
result.teta=teta;
result.ab=ab;
