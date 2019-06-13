%20180523
%ED full spectrum by diagonalazation the Floquet Hamiltonian in each sector
%using symmetries to speed up
% get the matrix element of X
tic
N_atom=12;
epslist=0.0002:0.0002:0.001;
T=0.005;
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

D=cell(2,length(epslist));

H_int=-Hamiltonian(N_atom,bc,cplist,'dipolar',1);

piz=PiPulse(N_atom,3);
seclist=[find(diag(piz)==1),find(diag(piz)==-1)];

[xx,yy]=meshgrid(epslist,1:2^(N_atom-1));
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
xelelist=zeros(size(epslist));
ndglist=xelelist;
for count=1:length(epslist)
    xele=0;
    ndg=0; % numer of degenerate pairs
    eps=epslist(count);
    Hp=Sigma_z*eps/2;
    for sym1=1:length(pt)/2
        sym2=sym1+length(pt)/2;
        
        Dtemp=[];
        U1=pt{sym1};
        
        if isempty(U1)
            continue
        end
        U2=pt{sym2};
        
        Hint1=U1'*H_int*U1;
        Hint2=U2'*H_int*U2;
        Hp1=U1'*Hp*U1;
        Hp2=U2'*Hp*U2;
        
        X12=U1'*Sigma_x*U2;
        [V1,D1]=eig(full(Hint1));
        [V2,D2]=eig(full(Hint2));
        D1=sparse(diag(D1));
        D2=sparse(diag(D2));
        %         [Vp1,Dp1]=eig(full(Hp1));
        %         [Vp2,Dp2]=eig(full(Hp2));
        Uf1=V1*diag(exp(1i*D1*T))*V1'*diag(exp(1i*sparse(diag(Hp1))));
        Uf2=V2*diag(exp(1i*D2*T))*V2'*diag(exp(1i*sparse(diag(Hp2))));
        [Vf1,Df1]=eig(full(Uf1));
        [Vf2,Df2]=eig(full(Uf2));
        Df1=1i*log(Df1);
        Df2=1i*log(Df2);
        Df1=diag(Df1);
        Df2=diag(Df2);
        
        for p1=1:length(Df1)
            pp1=Df1(p1);
            for p2=1:length(Df2)
                pp2=Df2(p2);
                if abs(pp1-pp2)<1e-4
                    xele=xele+abs(Vf1(:,p1)'*X12*Vf2(:,p2))^2;
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

figure(6)
plot(10*epslist,xelelist)
xlabel('\epsilon')
ylabel('ungapped X elements')
ppStyle(20,2)
hold on
figure(7)
plot(epslist,ndglist)
xlabel('\epsilon')
ylabel('percentage of no-gap states')
ppStyle(20,2)
hold on
% filetosave='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/ED_spec_sym';
% for p = fix(clock)
%     filetosave=[filetosave,'_',num2str(p)];
% end
% savefig([filetosave,'.fig'])
% save([filetosave,'.mat'],'N_atom','glist','bc','cplist','D','gap')