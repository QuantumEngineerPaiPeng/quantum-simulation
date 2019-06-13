%20180106
%Compute ZZ OTOC via diagonalizing in fermion picture
%H=-gZ-JXX
%Note here H is not the Hamiltonian, it's defined by the propagator in the
%fermion picture: U=e^(iHt). H here is twice the Hamiltonian.
N_spin=10;
J=1;
%t=0.16*12;

rho0=diag([ones(1,N_spin),-ones(1,N_spin)]);
rho0=rho0/(2*N_spin)^0.25;

tlist=10:0.2:40;
glist=0:0.1:3;

OTOC=zeros(length(glist),length(tlist));
TPC=zeros(length(glist),length(tlist));

row=1;
for g=glist
    H11=diag(ones(1,N_spin-1),1);
    %H11(1,N_spin)=-1;
    H12=H11-transpose(H11);
    %H12(1,N_spin)=1;
    %H12(N_spin,1)=-1;
    H11=H11+transpose(H11);
    H=-J/2*[H11,H12; % XX
        -H12,-H11];
%     H=-J*[H11,zeros(N_spin); % XX+YY
%         zeros(N_spin),H11];
    H=H+g*diag([-ones(1,N_spin),ones(1,N_spin)]);
    
    [V,D]=eig(H);
    rhoV=transpose(V)*rho0*V;
    
    col=1;
    for t=tlist
        rhot=V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*transpose(V);
        OTOC(row,col)=-trace((rhot*rho0-rho0*rhot)^2);
        TPC(row,col)=trace(rhot*rho0);
        col=col+1;
    end
    row=row+1;
end
% [maxC,temp]=max(OTOC,[],1);
% maxg=glist(temp);
% figure(1)
% hold on
% plot(tlist,maxg)
% figure(2)
% hold on
% plot(tlist,maxC)
figure(2)
hold on
plot(glist/2,2*mean(real(OTOC),2))
% figure(3)
% hold on
% plot(glist/2,2*(real(TPC))
