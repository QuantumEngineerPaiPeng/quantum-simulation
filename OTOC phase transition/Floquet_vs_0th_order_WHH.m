%20180322
%simulation of the WAHUHA sequence considering the finite width of the
%pulse to compare the Floquet Hamiltonian and propagator with the 0th
%order case
J=8.18e-3;
p1=1;
Tlist=logspace(1,3,20);
N_atom=8;

NT=100;
u=-0.2;
v=0.2;
w=0;

bc='p';% 'p' for periodic boundary condition; 'o' for open boundary condition
% cplist=[1,1/8,1/27];
cplist=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
for pp=1:length(cplist)
    for p=1:N_atom-1
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+J*cplist(pp)*...
            (Sigma_i{3}{p}*Sigma_i{3}{q}-...
            (Sigma_i{1}{p}*Sigma_i{1}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
    end
end

Hdip=zeros(2^N_atom);
for pp=1:length(cplist)
    for p=1:N_atom-1
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        Hdip=Hdip+J*cplist(pp)*...
            (Sigma_i{2}{p}*Sigma_i{2}{q}-...
            (Sigma_i{1}{p}*Sigma_i{1}{q}+Sigma_i{3}{p}*Sigma_i{3}{q})/2);
    end
end

px=expm(1i*(H_int-(pi+0.1)/(4*p1)*Sigma_x)*p1);% here are U^\dag
py=expm(1i*(H_int-(pi+0.1)/(4*p1)*Sigma_y)*p1);
pxb=expm(1i*(H_int+(pi+0.1)/(4*p1)*Sigma_x)*p1);
pyb=expm(1i*(H_int+(pi+0.1)/(4*p1)*Sigma_y)*p1);

[V,D]=eig(H_int);

Hfid=[];
Ufid=[];

for cycle=Tlist
tau=cycle/6;

t1=tau-p1/2;
t2=tau-p1;

Uf=V*diag(exp(diag(1i*D*t1)))*V'*px*V*diag(exp(diag(1i*D*t2)))*V'*py*...
    V*diag(exp(diag(1i*2*D*t1)))*V'*pyb*V*diag(exp(diag(1i*D*t2)))*V'*...
    pxb*V*diag(exp(diag(1i*D*t1)))*V';

Hf=logm(Uf)/1i/cycle;

Hfid=[Hfid;trace(Hf*Hf')];
Ufid=[Ufid;trace(Uf^NT)];
end

figure(1)
plot(Tlist,Hfid)
figure(2)
hold on
plot(Tlist,abs(Ufid))