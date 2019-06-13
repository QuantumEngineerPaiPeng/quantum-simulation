% 20180821
% ED full spectrum by diagonalazation and get the entanglement entropy
% Created by Pai Peng
tic
N_atom=10;
% glist=logspace(-2,0,5);
% glist=linspace(0.01,10,100);
g=3;
bc='p';
model='dq';
cplist=(1:5).^-3;
% cplist=1;

H_int=OperatorClass(N_atom,model,-1,bc,cplist,1);
Hz=OperatorClass(N_atom,'z',1);

% [xx,yy]=meshgrid(glist,1:2^(N_atom-1));
EE=zeros(1,2^N_atom);
H=H_int+g*Hz;
H.diagonalize()
E=H.eigsys{1}.D;
for count=1:2^N_atom
    psi=H.eigsys{1}.V(:,count);
    rho=OperatorClass(N_atom);
    rho.matrix={psi*psi'};
    EE(count)=BipartiteEE(rho,5);
end
[E,ind]=sort(E);
Olap=Olap(ind);

% figure
% % scatter(E,EE,'*')
% scatter(1:2^N_atom,EE,'*')
% % set(gca,'yscale','log')
% xlabel('E_n')
% ylabel('Entanglement entropy')
% % ylim([1e-10,1])
% ppStyle(20,1,5)

toc

% filetosave='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ED_spec_sym';
% for p = fix(clock)
%     filetosave=[filetosave,'_',num2str(p)];
% end
% savefig([filetosave,'.fig'])
% save([filetosave,'.mat'],'N_atom','glist','bc','cplist','D','gap')