%20180406
%Compute tr([U1,U2]^2) with U1, U2 being t-Ising model with field strength
%g1 and g2
%Note here H is not the Hamiltonian, it's defined by the propagator in the
%fermion picture: U=e^(iHt). H here is twice the Hamiltonian.
tic
N_spin=50;
J=1;
%t=0.16*12;


tlist=40:0.1:160;
glist=logspace(-1,1,21);

[X,Y]=meshgrid(glist);
OTOC=zeros(size(X));

H11=diag(ones(1,N_spin-1),1);
%H11(1,N_spin)=-1;
H12=H11-transpose(H11);
%H12(1,N_spin)=1;
%H12(N_spin,1)=-1;
H11=H11+transpose(H11);
H=-J/2*[H11,H12; % XX
    -H12,-H11];

row=1;
for g1=glist
    col=1;
    for g2=glist
        if g2<=g1
            OTOC(row,col)=OTOC(col,row);
            col=col+1;
            continue
        end
        H1=H+g1*diag([-ones(1,N_spin),ones(1,N_spin)]);
        H2=H+g2*diag([-ones(1,N_spin),ones(1,N_spin)]);
        [V1,D1]=eig(H1);
        [V2,D2]=eig(H2);
        for t=tlist
            U1=V1*diag(exp(diag(-1i*D1*t)))*V1';
            U2=V2*diag(exp(diag(-1i*D2*t)))*V2';
            
            OTOC(row,col)=OTOC(row,col)+trace((U1*U2-U2*U1)*(U1*U2-U2*U1)')/length(tlist);
        end
        col=col+1;
    end
    row=row+1;
end


figure(3)
h=pcolor(X,Y,real(OTOC));
set(h, 'EdgeColor', 'none');
colorbar
toc