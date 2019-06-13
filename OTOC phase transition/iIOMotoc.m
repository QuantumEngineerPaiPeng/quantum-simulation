%20180214
%Compute OTOC for iIOM
N_spin=50;
J=1;
g_list=0.:0.5:1.5;
tlist=12:0.1:40;
%t=0.16*12;
otoclist=[];

for g=g_list
    H11=diag(ones(1,N_spin-1),1);
    %H11(1,N_spin)=-1;
    H12=H11-transpose(H11);
    %H12(1,N_spin)=1;
    %H12(N_spin,1)=-1;
    H11=H11+transpose(H11);
    H=-J/2*[H11,H12;
        -H12,-H11];
    
    H=H+g*diag([-ones(1,N_spin),ones(1,N_spin)]);
    
    [V,D]=eig(H);
    
    Z0=diag([ones(1,N_spin),-ones(1,N_spin)]);
    %rho0=rho0/(2*N_spin)^0.25;
    otoc=0;
    for p=1:N_spin
        rho0=zeros(size(Z0));
        rho0(p,p)=1;
        rho0(N_spin+p,N_spin+p)=-1;
        rho0=rho0/sqrt(2*N_spin);
        
        rhoV=transpose(V)*rho0*V;
        
        rhobar=zeros(size(rho0));
        for t=tlist
            rhobar=rhobar+V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*transpose(V);
        end
        rhobar=rhobar/(length(tlist)+1);
        otoc=otoc-trace((rhobar*Z0-Z0*rhobar)^2)/2*4;%/trace(rhobar*rhobar);
    end
    otoc=otoc/N_spin;
    otoclist=[otoclist,otoc];
end

figure(1)
hold on
plot(g_list,otoclist,'LineWidth',2)