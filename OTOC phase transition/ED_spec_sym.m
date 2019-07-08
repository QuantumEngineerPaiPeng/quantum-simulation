%20180422
%ED full spectrum by diagonalazation in each sector
%using symmetries to speed up
tic
N_atom=10;
% glist=logspace(-2,0,5);
glist=linspace(0.01,10,100);
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
cplist=[1,1/8,1/27];
% cplist=1;

D=cell(2,length(glist));

H_int=-Hamiltonian(N_atom,bc,cplist,'Ising',1);

piz=PiPulse(N_atom,3);
seclist=[find(diag(piz)==1),find(diag(piz)==-1)];

% [xx,yy]=meshgrid(glist,1:2^(N_atom-1));
toc
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
if exist(filename,'file')
    load(filename)
else
    pt=AllkpzProject(N_atom);
    save(filename,'pt')
end
toc
count=1;
parfor count=1:length(glist)
    g=glist(count);
    H=H_int+g*Sigma_z;
    for sec=1:2
%         Hs=H(seclist(:,sec),seclist(:,sec));
        Dtemp=[];
        symlist=(1:length(pt)/2)+(sec-1)*length(pt)/2;
        for sym=symlist
            U=pt{sym};
            if isempty(U)
                continue
            end

            Hsym=U'*H*U;
            Dtemp=[Dtemp;eig(full(Hsym))];
        end
        D{sec,count}=[D{sec,count},sort(Dtemp)];
    end
end
D=cell2mat(D);
D1=D(1:end/2,:);
D2=D(end/2+1:end,:);
gap=abs(D1-D2);
Dmax=(0.5*(D1(end,:)+D2(end,:)))/N_atom;
Dmin=(0.5*(D1(1,:)+D2(1,:)))/N_atom;
gap=gap./Dmax;
Dmean=real((0.5*(D1+D2)))/N_atom;
Drel=(Dmean-Dmin)./(Dmax-Dmin);
[xx,yy]=meshgrid(glist,1:2^(N_atom-1));
toc

% for p=1:length(glist)
%     figure
%     histfit(gap(:,p),2^N_atom/2,'exponential')
% end
figure
pcolor(xx,yy,abs(gap))
colorbar
shading interp
% filetosave='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ED_spec_sym';
% for p = fix(clock)
%     filetosave=[filetosave,'_',num2str(p)];
% end
% savefig([filetosave,'.fig'])
% save([filetosave,'.mat'],'N_atom','glist','bc','cplist','D','gap')