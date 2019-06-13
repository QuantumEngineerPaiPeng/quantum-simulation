%20180422
%Eigenvalue statics for dipolar model
tic
N_atom=14;
T=0.5;
epslist=0.0:0.03:0.3;
bc='p';
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
%RandomM=rand(N_rand,3,N_atom);
cplist=[1,1/8,1/27];
r=[];
toc
H_int=Hamiltonian(N_atom,bc,cplist,'dipolar',1);

parity=ParityOperator(N_atom);
piz=PiPulse(N_atom,3);
k0=kProject(N_atom,0);
translation=TranslationOperator(N_atom,1);
toc
H_intk=full(k0'*H_int*k0);
Z_k=full(k0'*Sigma_z*k0);
U_int=expm(-1i*H_intk*T);

row=1;
for eps=epslist
    U_p=expm(-1i*Z_k*(1+eps)*pi/2);
    U=U_p*U_int;
    [V,D]=eig(U);
    D=real(1i*log(D));
%     [V,D]=eig(full(k0'*(H_int+pi*Sigma_z)*k0));
    D=diag(D);
    PP=diag(V'*k0'*parity*k0*V);
    ZZ=diag(V'*k0'*piz*k0*V);
    sec=[];
    for p=1:length(PP)
        if and(abs(PP(p)-1)<1e-10,abs(ZZ(p)-1)<1e-10)
            sec=[sec,p];
        end
    end
    Dd=sort(D(sec));
%     Dd=sort(D);
    delta=Dd(2:end)-Dd(1:end-1);
    for j=1:(length(delta)-1)
        r(row,j)=min([delta(j),delta(j+1)])/max([delta(j),delta(j+1)]);
    end
    row=row+1;
end
toc
r_ave=real(mean(r,2));
x=0:0.01:1;
% figure(1)
% histogram(r,20)
% hold on
% plot(x,2*(1+x).^(-2)*length(r)/20)
% plot(x,27/4*(x+x.^2).*(1+x+x.^2).^(-2.5)*length(r)/20)
figure(2)
plot(epslist,r_ave)