%20180523
%ED full spectrum by diagonalazation in each sector
%using symmetries to speed up
% get the matrix element of X
tic
N_atom=15;
glist=0.05:0.05:0.5;
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

D=cell(2,length(glist));

H_int=-Hamiltonian(N_atom,bc,cplist,'dipolar',1);

piz=PiPulse(N_atom,3);
seclist=[find(diag(piz)==1),find(diag(piz)==-1)];

[xx,yy]=meshgrid(glist,1:2^(N_atom-1));
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
xelelist=zeros(size(glist));
ndglist=xelelist;
parfor count=1:length(glist)
    xele=0;
    ndg=0; % numer of degenerate pairs
    g=glist(count);
    H=H_int+g*Sigma_z;
    for sym1=1:length(pt)/2
        sym2=sym1+length(pt)/2;
        
        
        Dtemp=[];
        U1=pt{sym1};
        
        if isempty(U1)
            continue
        end
        U2=pt{sym2};
        
        Hsym1=U1'*H*U1;
        Hsym2=U2'*H*U2;
        X12=U1'*Sigma_x*U2;
        [V1,D1]=eig(full(Hsym1));
        [V2,D2]=eig(full(Hsym2));
        D1=diag(D1);
        D2=diag(D2);
        for p1=1:length(D1)
            pp1=D1(p1);
            for p2=1:length(D2)
                pp2=D2(p2);
                if abs(pp1-pp2)<1e-4
                    xele=xele+abs(V1(:,p1)'*X12*V2(:,p2))^2;
%                     xele=xele+V1(:,p1)'*X12*V2(:,p2);
                    ndg=ndg+1;
                end
            end
        end
    end
    xelelist(count)=xele;
    ndglist(count)=ndg;
end
xelelist=xelelist/(N_atom*2^N_atom);
ndglist=ndglist/(2^N_atom/2);
toc

figure(1)
plot(glist,xelelist)
xlabel('g')
ylabel('gapped X elements')
ppStyle(20,2)
hold on
figure(2)
plot(glist,ndglist)
xlabel('g')
ylabel('percentage of no-gap states')
ppStyle(20,2)
hold on
% filetosave='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ED_spec_sym';
% for p = fix(clock)
%     filetosave=[filetosave,'_',num2str(p)];
% end
% savefig([filetosave,'.fig'])
% save([filetosave,'.mat'],'N_atom','glist','bc','cplist','D','gap')