function y=Ising_gsMz_integr(q,g)
wq=sqrt(g.^2+2*g.*cos(q)+1)./g;
y=1/pi*(1./wq+1./g.*cos(q)./wq);
end