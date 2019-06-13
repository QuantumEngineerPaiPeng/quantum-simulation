% 20180811
% ED full spectrum by diagonalazation in each sector and get the overlap
% of the eigenvector with a state |\psi>
% Using symmetries to speed up. Note that rho does not need to obey the
% symmetry.
% Created by Pai Peng
tic
N_atom=12;
% glist=logspace(-2,0,5);
% glist=linspace(0.01,10,100);
g=10;
bc='p';
model='dipolar';
sym='kpz';
cplist=(1:5).^-3;
% cplist=1;

D=cell(2,length(glist));

H_int=OperatorClass(N_atom,model,-1,bc,cplist,3);
H_int.symmetrize(sym);
Hz=OperatorClass(N_atom,'z',1);
Hz.symmetrize(sym);
dir1=[90,0];
dir2=[90,0];
rho=StateN({dir1,dir2,dir1,dir2,dir1,dir2,dir1,dir2,dir1,dir2,dir1,dir2});
% rho.matrix={zeros(2^N_atom,2^N_atom)};
% rho.matrix{1}(1,1)=1;
rho.symmetrize(sym);

% [xx,yy]=meshgrid(glist,1:2^(N_atom-1));
E=[];
Olap=[];

H=H_int+g*Hz;
H.diagonalize()

for p=1:length(H.matrix)
    E=[E,H.eigsys{p}.D.'];
    for q=1:length(H.eigsys{p}.D)
        Olap=[Olap,H.eigsys{p}.V(:,q)'*rho.matrix{p}*H.eigsys{p}.V(:,q)];
    end
end
[E,ind]=sort(E);
Olap=Olap(ind);


figure
scatter(E,Olap,'*')
% scatter(1:2^N_atom,Olap,'*')
set(gca,'yscale','log')
xlabel('E_n')
ylabel('|<X|\psi_n>|^2')
ylim([1e-10,1])
ppStyle(20,1,5)

toc

% filetosave='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ED_spec_sym';
% for p = fix(clock)
%     filetosave=[filetosave,'_',num2str(p)];
% end
% savefig([filetosave,'.fig'])
% save([filetosave,'.mat'],'N_atom','glist','bc','cplist','D','gap')