%20180427
%OTOC for different t and g
%using symmetries to speed up
% function version for cluster
function OTOC_t_g_sym_func(N_atom,glist,tlist,bc,cplist,model)
tic
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


H_int=-Hamiltonian(N_atom,bc,cplist,model,1);

piz=PiPulse(N_atom,3);

toc
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',N_atom);
if exist(filename,'file')
    load(filename)
else
    pt=AllkpzProject(N_atom);
    save(filename,'pt')
end
toc

rho_0=Sigma_z;
V0=Sigma_x;

C=zeros(length(pt)/2,length(glist),length(tlist));
for sym=1:length(pt)/2
    Up=pt{sym};
    Um=pt{sym+length(pt)/2};
    if or(isempty(Up),isempty(Um))
        continue
    end
    rhop=Up'*rho_0*Up;
    rhom=Um'*rho_0*Um;
    V0pm=Up'*V0*Um;
    parfor count=1:length(glist)
        g=glist(count);
        H=H_int+g*Sigma_z;
        
        Hp=Up'*H*Up;
        Hm=Um'*H*Um;
        [Vp,Dp]=eig(full(Hp));
        [Vm,Dm]=eig(full(Hm));
        Dp=sparse(diag(Dp));
        Dm=sparse(diag(Dm));
        Ctemp=[];
        for t=tlist
            rhopt=Vp*diag(exp(1i*Dp*t))*Vp'*rhop*Vp*diag(exp(-1i*Dp*t))*Vp';
            rhomt=Vm*diag(exp(1i*Dm*t))*Vm'*rhom*Vm*diag(exp(-1i*Dm*t))*Vm';
            cm=2*sum(sum(abs(rhopt*V0pm-V0pm*rhomt).^2));
            Ctemp=[Ctemp,cm];
        end
        C(sym,count,:)=Ctemp;
    end
    save('temporary.mat','C','sym')
end

OTOC=squeeze(sum(C,1))/(4*N_atom*2^N_atom);
toc
[maxC,temp]=max(OTOC,[],1);
maxg=glist(temp);

pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
for p = fix(clock)
    filetosave=[filetosave,num2str(p),'_'];
end
save([filetosave,'.mat'],'N_atom','glist','tlist','bc','cplist','OTOC','model')
figure(1)
plot(tlist,maxg)
savefig([filetosave,'fig1','.fig'])
figure(2)
plot(tlist,maxC)
savefig([filetosave,'fig2','.fig'])
figure(3)
plot(glist,OTOC)
savefig([filetosave,'fig3','.fig'])

end