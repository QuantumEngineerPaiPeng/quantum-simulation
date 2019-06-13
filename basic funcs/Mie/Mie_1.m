function result = Mie_1(m, x, xmax)
% Rayleigh Approximation of Mie Efficiencies for given 
% complex refractive-index ratio m=m'+im" 
% and size parameter x=k0*a, k0= wave number in ambient medium,  
% a=sphere radius using the low-frequency approximation of the 
% complex Mie Coefficients an and bn for n=1, p. 131 
% in Bohren and Huffman (1983) BEWI:TDD122.
% Result is m', m", x, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% asymmetry parameter asy=<costeta> and qratio=qb/qsca.
% Uses the function "mieab_1" for the Mie Coefficient a1.
% C. Mätzler, May 2002.

if x==0,
    result=[0 0 0 0 0 1.5];
elseif x>xmax,
    result=[NaN, NaN, NaN, NaN, NaN, NaN];
else,
    nmax=1; n=1; cn=2*n+1; x2=x*x;
    f=mieab_1(m,x);
    anp=real(f); anpp=imag(f);
    dn=cn.*anp;
    qabs=2*dn/x2;
    en=cn.*(anp.*anp+anpp.*anpp);
    qsca=2*en/x2;   
    qext=qabs+qsca;
    q=-cn.*f;
    qb=q*q'/x2;
    qratio=qb/qsca;
    result=[qext qsca qabs qb 0 qratio];
end;