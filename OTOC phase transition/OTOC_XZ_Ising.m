%20180121
%Compute CXZ and CYZ OTOC for intergrable model
N_spin=100;
J=1;
tlist=0:0.01:0.5;
glist=.0:0.1:3;
%t=0.16*12;
cnst=0.5;
H11=diag(ones(1,N_spin-1),1);
% H11(1,N_spin)=-1;
H12=H11-transpose(H11);
% H12(1,N_spin)=1;
% H12(N_spin,1)=-1;
H11=H11+transpose(H11);
H_int=-J/2*[H11,H12; % XX
    -H12,-H11];
% H_int=-J*[H11,cnst*H12; % XX+YY+cnst(XX-YY)
%     -cnst*H12,-H11];
rho0=diag([ones(1,N_spin),-ones(1,N_spin)]);
rho0=rho0/(2*N_spin)^0.5;

CZX=zeros(length(glist),length(tlist));
CZY=zeros(length(glist),length(tlist));

row=1;
for g=glist
    col=1;
    H=H_int+g*diag([-ones(1,N_spin),ones(1,N_spin)]);
    [V,D]=eig(H);
    rhoV=transpose(V)*rho0*V;
    for t=tlist
        rhot=V*diag(exp(diag(1i*D*t)))*rhoV*diag(exp(diag(-1i*D*t)))*transpose(V);
        rho11=rhot(1:N_spin,1:N_spin);
        rho12=rhot(1:N_spin,N_spin+1:2*N_spin);
        
        CZX(row,col)=sum(diag(rho11).^2);
        CZY(row,col)=sum(diag(rho11).^2);
        
        M1=(real(rho11+rho12)).^2;
        M2=(real(rho11-rho12)).^2;
        M3=(imag(rho11)).^2+(imag(rho12)).^2;
        
        for k=1:N_spin-1
            CZX(row,col)=CZX(row,col)+(k+1)*sum(diag(M2,k))+(k-1)*sum(diag(M1,k))+2*k*sum(diag(M3,k));
            CZY(row,col)=CZY(row,col)+(k+1)*sum(diag(M1,k))+(k-1)*sum(diag(M2,k))+2*k*sum(diag(M3,k));
        end
        col=col+1;
    end
    row=row+1;
end

[maxC,temp]=max(CZX,[],1);
maxg=glist(temp);
figure(1)
hold on
plot(tlist*2,maxg/2)
figure(2)
hold on
plot(tlist*2,maxC)
figure(3)
hold on
plot(glist/2,CZX)
% figure(2)
% hold on
% dotoc=CZX(:,2:size(CZY,2))-CZX(:,1:size(CZY,2)-1);
% plot(tlist(1:length(tlist)-1),dotoc/(tlist(2)-tlist(1)))