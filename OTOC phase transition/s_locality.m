% 20180817
% Using the saved prethermal Hamiltonian to calculate locality of s
% as a function of epsilon
tic
N_atom=8;
fname=fullfile(...
    '~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/prethermalH_func'...
    ,['s',num2str(N_atom)]);
load(fname)
% order=length(hpre)-1;

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
% epslist=logspace(-1,1,10);
epslist=0.5:0.5:3;
Cl=zeros(length(epslist),N_atom+1);

pfolder='~/Dropbox (MIT)/grad/research/codes/OTOC phase transition/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];




for count=1:length(epslist)
    eps=epslist(count); % eps is the prefactor of interaction in spin picture
    so=speye(length(s{1}));
    
    for pp=1:length(s)
        so=so+eps^(pp)*s{pp};
    end
    Cl(count,:)=CLWeight(so);
    count
end

%     save([filetosave,'r_',num2str(N_atom)]);
toc


figure(1)
semilogy(0:N_atom,Cl)
xlabel('order')
ylabel('r')
ppStyle(20,2)