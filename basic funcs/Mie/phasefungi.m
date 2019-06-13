function result = phasefungi(m,x)

% Legendre coefficients gi of phase function p(costeta) 
% in Mie Theory (e.g. gl in
% Meador and Weaver, 1980 J. Atm. Sci. 37, pp. 630-643)

% Input: 
% m: refractive index, x: size parameter
% nj: number of gi values (use about nj=2.5*(x+1))
% nsteps: number of angular values between 0 and pi
% (take about nsteps = 20*nj)
% C. Mätzler, July 2003

nj=ceil((x+1.5)*2.3)
nsteps=25*nj;
dteta=pi/nsteps;
m1=real(m); m2=imag(m);
nx=(1:nsteps);
teta=(nx-0.7)*dteta; u=cos(teta); s=sin(teta);
    for j = 1:nsteps, 
        a(:,j)=Mie_S12(m,x,u(j));
        SL(j)= real(a(1,j)'*a(1,j));
        SR(j)= real(a(2,j)'*a(2,j));
        A(j)=(SL(j)+SR(j)).*s(j);
    end;
Q=mie(m,x);
Qext=Q(1); Qsca=Q(2); asy=Q(5); w=Qsca/Qext;
AA=dteta*A./Qsca/x.^2;

for jj=1:nj,
    x=legendre(jj-1,u);
    x0=x(1,:);
    gi(jj)=x0*AA';
end;
result=gi;