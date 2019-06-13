%20170312
%analytically compute the ground state energy of Ising model with
%transverse field
% glist=0.:0.01:2;
E0=[];
E0p=[];
for g=glist
    E0=[E0,integral(@(k) Ising_gsMz_integr(k,g),0,pi)];
    E0p=[E0p,integral(@(k) Ising_gsMz_integr(k,g+1e-5),0,pi)];
end

figure(2)
hold on
plot(glist/2,E0)
figure(3)
hold on
plot(glist/2,(E0p-E0)*1e5)
       