%20180308
% compute the ground state using exact diagonalization
% Using H/phi to iteratively compute the ground state
% Using only matrix multiplication
% Use symmetry to speed up
tic
N_atom=14;
N_rand=1;
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
cplist=[1,1/8,1/27];
% cplist=1;
model='dipolar';
bc='p';
glist=0.2:0.2:3;
tlist=0.5:0.2:4;
clear('Sigma_i')
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
if exist(filename,'file')
    load(filename)
else
    pt=AllkpzProject(N_atom);
    save(filename,'pt')
end
sym={pt{length(pt)/2+2},pt{2}};
H_int=-Hamiltonian(N_atom,bc,cplist,model,1);
V1=Sigma_z;
V2=Sigma_x;

Y2=[];
energy=zeros(length(glist),2);
Csec=zeros(length(glist),length(tlist),2);

v1={sym{1}'*V1*sym{1},sym{2}'*V1*sym{2}};
v2={sym{1}'*V2*sym{2},(sym{1}'*V2*sym{2})'};

for row=1:length(glist)
    energytemp=zeros(1,2);
    Csectemp=zeros(length(tlist),2);
    H=H_int+glist(row)*Sigma_z;
    h={sym{1}'*H*sym{1},sym{2}'*H*sym{2}};
    V={};
    D={};
    gs={};
    
    for sec=1:2
        proj=sym{sec};
        dim=size(proj,2);
        [Vtemp,Dtemp]=eig(full(h{sec}));
        Dtemp=sparse(diag(Dtemp));
        [E,p]=min(Dtemp);
        energytemp(sec)=E;
        gstemp=Vtemp(:,p);
        V{end+1}=Vtemp;
        D{end+1}=Dtemp;
        gs{end+1}=gstemp;
        
    end
    
    energy(row,:)=energytemp;
    
    for col=1:length(tlist)
        t=tlist(col);
        v1t={};
        for sec=1:2
            v1t{end+1}=V{sec}*diag(exp(1i*D{sec}*t))*V{sec}'...
                *v1{sec}*...
                V{sec}*diag(exp(-1i*D{sec}*t))*V{sec}';
        end
        commutator={(v1t{1}*v2{1}-v2{1}*v1t{2})*(v2{2}*v1t{1}-v1t{2}*v2{2}),...
            (v2{2}*v1t{1}-v1t{2}*v2{2})*(v1t{1}*v2{1}-v2{1}*v1t{2})};
        for sec=1:2
            Csectemp(col,sec)=gs{sec}'*commutator{sec}*gs{sec};
        end
    end
    Csec(row,:,:)=Csectemp;
end
Csec=Csec/(2^N_atom);
OTOC=squeeze(sum(Csec,3));
[maxC,temp]=max(OTOC,[],1);
maxg=glist(temp);
figure(1)
plot(tlist,maxg)
figure(2)
plot(tlist,maxC)
figure(3)
plot(glist,OTOC)
toc