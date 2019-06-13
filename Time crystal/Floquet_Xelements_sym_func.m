%20180523
%ED full spectrum by diagonalazation the Floquet Hamiltonian in each sector
%using symmetries to speed up
% get the matrix element of X
function Floquet_Xelements_sym_func(N_atom,epslist,T,bc,cplist,model,thr)
tic
% N_atom=14;
% epslist=0.02:0.02:0.1;
% T=0.5;
% bc='p';
Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
% cplist=[1,1/8,1/27];



H_int=-Hamiltonian(N_atom,bc,cplist,model,1);

toc
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
if exist(filename,'file')
    load(filename)
else
    pt=AllkpzProject(N_atom);
    save(filename,'pt')
end
toc

pfolder='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
for p = fix(clock)
    filetosave=[filetosave,num2str(p),'_'];
end

xelelist=zeros(size(epslist));
ndglist=xelelist;
parfor count=1:length(epslist)
    xele=0;
    ndg=0; % numer of degenerate pairs
    eps=epslist(count);
    Hp=Sigma_z*eps/2;
    for sym1=1:length(pt)/2
        sym2=sym1+length(pt)/2;
        
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
                if abs(pp1-pp2)<thr
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

figure(1)
plot(epslist,xelelist)
xlabel('\epsilon')
ylabel('gapped X elements')
ppStyle(20,2)
hold on
figure(2)
plot(epslist,ndglist)
xlabel('\epsilon')
ylabel('percentage of no-gap states')
ppStyle(20,2)
hold on

save([filetosave,'.mat'],'N_atom','epslist','bc','cplist','model','thr','xelelist','ndglist','T')
writelog(pfolder,{'N_atom','epslist','bc','cplist','model','thr','xelelist','ndglist','T'},...
    {N_atom,epslist,bc,cplist,model,thr,xelelist,ndglist,T})