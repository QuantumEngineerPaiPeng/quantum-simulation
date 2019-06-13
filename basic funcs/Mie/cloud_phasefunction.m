function result = cloud_phasefunction(lamda, rc, mu, nsteps)

% Computation of Phase Function pm (unpolarised)
% of clouds with normalisation to 'one' when integrated 
% over all directions/(4*pi), see Chandrasekhar 1960
% Eq. (28), 
% Input: complex refractive index m=m'+im", 
% mode radius rc, wavelength lamda, both in micron, 
% mu=cos(scattering angle), and nsteps=number of drop sizes  
% Output: pm,  
% asym=asymmetry parameter, and wm=single-scattering albedo
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, July 2003

dr0=0.5*rc;
dr=1.5*rc/nsteps;
mv=1; alfa=6; gam=1;   % Cloud parameters according to Ulaby et al.1981
m=nwater(lamda); 
nx=(1:nsteps)';
r=(nx-1)*dr+dr0;
rm=1e-6*r; rcm=1e-6*rc; drm=dr*1e-6;   % transformations to m
xj=2*pi*r./lamda;       % size parameter
sigmag=pi*rm.*rm;      % geometric cross section in m2
NMP=modgammaulaby(rm,rcm,mv,alfa,gam); 
for j=1:nsteps,        % computations for increasing drop size     
    x=xj(j);
    nmax=round(2+x+4*x^(1/3));
    ab=mie_ab(m,x);
    an=ab(1,:);
    bn=ab(2,:);
    pt=mie_pt(mu,nmax);
    pin =pt(1,:);
    tin =pt(2,:);
    n=(1:nmax);
    n2=(2*n+1)./(n.*(n+1));
    pin=n2.*pin;
    tin=n2.*tin;
    S1=(an*pin'+bn*tin');
    S2=(an*tin'+bn*pin');
    Q=mie(m,x);
    Qext(j)=NMP(j)*Q(1); Qsca=Q(2); asy(j)=NMP(j)*Q(5); 
    w(j)=NMP(j)*Qsca/Q(1);
    p(j)=2*NMP(j)*(S1'*S1+S2'*S2)/(Qsca*x^2);  
end;
sNMP=sum(NMP);
pm=sum(p)/sNMP;
asym=sum(asy)/sNMP;
wm=sum(w)/sNMP;
Qextm=sum(Qext)/sNMP;
% Qsca to be exchanged by Qext above if normalisation to 
% single-scattering albedo, w, is required
result=[pm,wm,asym,Qextm];