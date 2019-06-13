%20180422
%ED full spectrum by diagonalazation in each sector
tic
N_atom=12;
glist=linspace(0.01,5,100);
% glist=2;
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
cplist=(1:3).^-3;

D1=[];
D2=[];

toc
H_int=Hamiltonian(N_atom,bc,cplist,'Ising',1);

piz=PiPulse(N_atom,3);
sec1=find(diag(piz)==1);
sec2=find(diag(piz)==-1);

toc
[xx,yy]=meshgrid(glist,1:2^(N_atom-1));
row=1;
for g=glist
    H=H_int+g*Sigma_z;
	H1=H(sec1,sec1);
    H2=H(sec2,sec2);
    D1=[D1,sort(eig(full(H1)))];
    D2=[D2,sort(eig(full(H2)))];
    row=row+1;
end
gap=abs(D1-D2);
gap=gap./(0.5*(D1(end,:)+D2(end,:)))*N_atom;
toc

figure
pcolor(xx,yy,abs(gap))
colorbar
shading interp