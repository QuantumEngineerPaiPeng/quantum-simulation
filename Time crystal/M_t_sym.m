%20180436
% Magnetizatino for different t
% using symmetries to speed up
tic
N_atom=12;
Nlist=logspace(0,7,50);
epslist=0.02:0.02:0.1;
epslist=2*epslist;
T=1;
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
clear('Sigma_i')
cplist=[1,1/8,1/27];

H_int=-Hamiltonian(N_atom,bc,cplist,'dipolar',1);

toc
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
if exist(filename,'file')
    load(filename)
else
    pt=AllkpzProject(N_atom);
    save(filename,'pt')
end
toc

rho_0=Sigma_x;
count=1;
M=zeros(length(pt)/2,length(epslist),length(Nlist));

for sym=1:length(pt)/2
    Up=pt{sym};
    Um=pt{length(pt)/2+sym};

    if or(isempty(Up),isempty(Um))
        continue
    end
    rho0pm=Up'*rho_0*Um;
    
    Hp=Up'*H_int*Up;
    Hm=Um'*H_int*Um;
    [Vp,Dp]=eig(full(Hp));
    [Vm,Dm]=eig(full(Hm));
    Dp=sparse(diag(Dp));
    Dm=sparse(diag(Dm));
    Zp=Up'*Sigma_z*Up;
    Zm=Um'*Sigma_z*Um;
    
    parfor p=1:length(epslist)
        eps=epslist(p);
        Ufp=Vp*diag(exp(1i*Dp*T))*Vp'*diag(exp(1i*diag(Zp)*eps/2));
        Ufm=Vm*diag(exp(1i*Dm*T))*Vm'*diag(exp(1i*diag(Zm)*eps/2));
        [Vfp,Dfp]=eig(full(Ufp));
        Dfp=sparse(diag(Dfp));
        [Vfm,Dfm]=eig(full(Ufm));
        Dfm=sparse(diag(Dfm));
        rho0V=Vfp'*rho0pm*Vfm;
        Ctemp=[];
        for N=Nlist

            rhotV=diag(Dfp.^N)*rho0V*diag(Dfm'.^N);
            cm=2*sum(sum(rhotV.*conj(rho0V)));
            Ctemp=[Ctemp,cm];
        end
        M(sym,p,:)=Ctemp;
    end
    
end

Mt=squeeze(sum(M,1))/(N_atom*2^N_atom);
toc

% pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
% foldername=mfilename;
% mkdir(pfolder,foldername)
% filetosave=[pfolder,foldername,'/'];
% for p = fix(clock)
%     filetosave=[filetosave,num2str(p),'_'];
% end
% save([filetosave,'.mat'],'N_atom','epslist','Nlist','bc','cplist','Mt')
figure(1)

semilogx(Nlist*T,Mt,'k')
hold on
% xlabel('# of Floquet periods')
xlabel('t')
ylabel('|M_z|')
