%20171031
%Compute integrable of motion of integrable model, diagonalizing in
%fermion picture
%H_Ising=-gZ-JXX
%H_Kitaev=-JXX for even and -gYY for odd
clear

model='Ising';
bc='open';
pl=true;

N_spin=16;
J=1;
% for pp=0.1:0.1:1.5
%     g=pp;
%t=0.16*12;
g=10;
if strcmp(model,'Kitaev')
    H11=zeros(N_spin,N_spin);
    H12=H11;
    for p=1:N_spin-1
        if mod(p,2)==1
            H12(p,p+1)=g;
            H11(p,p+1)=-g;
        else
            H12(p,p+1)=-J;
            H11(p,p+1)=-J;
        end
    end
    H11=H11+H11';
    H12=H12-H12';
    H=[H11,H12;H12',-H11];
end


if strcmp(model,'Ising')
    H11=diag(ones(1,N_spin-1),1);
    if strcmp(bc,'periodic')
        H11(1,N_spin)=-1;
    end
    H12=H11-transpose(H11);
    if strcmp(bc,'periodic')
        H12(1,N_spin)=1;
        H12(N_spin,1)=-1;
    end
    H11=H11+transpose(H11);
    H=-J/2*[H11,H12;
        -H12,-H11];
    H=H+g*diag([-ones(1,N_spin),ones(1,N_spin)]);
end
[V,D]=eig(H); %V'*H*V=D

%IOM=zeros(2*N_spin,2*N_spin,2*N_spin);

N_plot=ceil(sqrt(N_spin));

if pl
    figure(1)
end
%
A=0;
B=0;

for p=1:N_spin
    d1=zeros(2*N_spin,2*N_spin);
    d1(p,p)=-1;
    d1(2*N_spin+1-p,2*N_spin+1-p)=1;
    IOM=V*d1*V';
    %     PX=[zeros(N_spin,N_spin),eye(N_spin);
    %        eye(N_spin),zeros(N_spin,N_spin)];
    %     IOM=V(:,p)*V(:,p)'-PX*V(:,p)*V(:,p)'*PX;
    %rho11=IOM(1:N_spin,1:N_spin,p);
    rho11=IOM(1:N_spin,1:N_spin);
    %rho12=IOM(1:N_spin,(N_spin+1):(2*N_spin),p);
    rho12=IOM(1:N_spin,(N_spin+1):(2*N_spin));
    A=A+trace(rho11*rho11');
    B=B+trace(rho12*rho12');
        if pl
            subplot(N_plot,N_plot,p)
            image(abs(IOM).^2,'CDataMapping','scaled')
        end
end

% sz=[1,0;0,-1];
% SZ=kron(sz,eye(N_spin));
% U=zeros(2*N_spin,2*N_spin);
%
% for p=1:N_spin
%     U(p,p)=(-1)^p;
%     U(N_spin+p,N_spin+p)=(-1)^p;
% end

