% 20180802
% Perturbatively diagonalize a Heisenberg interaction with random field
% Created by Pai Peng
L=6;
w=10;
nM=6;
randZ=rand(L,1)*2-1;
% randZ=0.1:0.1:0.8;
H01=w*diag(randZ)+diag(ones(L-1,1),1)+diag(ones(L-1,1),-1);
[V,D]=eig(H01);
Z0=cell(1,L);

for p=1:L
    temp=zeros(L);
    temp(p,p)=1;
    tp=V'*temp*V;
    Z0{p}=Fermion2Pauli([tp,zeros(L);zeros(L),-tp]); % sigma_z^p in the diagonal basis
end

H0=Fermion2Pauli([D,zeros(L);zeros(L),-D]); % Noninteracting Hamiltonian in the diagonal basis
E0=full(diag(H0));

H1=sparse(2^L,2^L); % ZZ in the diagonal basis
for p=1:(L-1)
    H1=H1+Z0{p}*Z0{p+1};
end
H1=0.1*H1;

E=sort(eig(full(H0+H1)));
r=sum(abs(E-sort(eig(H0))));

s={};
Hpre=H0;
hpre={H0};
for p=1:nM
    h=expandSH(s,H0,p)+expandSH(s,H1,p-1);
    hpre{end+1}=diag(diag(h));
    Hpre=Hpre+hpre{end};
    temp=h-hpre{end};
    [row,col,val]=find(temp);
    for q=1:length(row)
        val(q)=-val(q)/(E0(col(q))-E0(row(q)));
    end
    s{end+1}=sparse(row,col,val,2^L,2^L);
    r(p+1)=sum(abs(E-sort(eig(full(Hpre)))));
end

figure(2)
semilogy(0:nM,r)
hold on