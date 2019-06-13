N_atomlist=16;
Nlist=1:20:2e3;
epslist=0.01:0.01:0.4;
T=0.5;
% cycle=36;
bc='p';
cplist=[1,2,3,4,5,6,7,8].^-2;
model='dipolar';
p2=2;
% % figure(1)
% % hold on
for N_atom=N_atomlist
M_t_sym_func(N_atom,Nlist,epslist,T,bc,cplist,model)
% PP9_func(N_atom,Nlist,epslist,cycle,bc,cplist,model,p2)
end
% pfolder='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
% foldername=mfilename;
% mkdir(pfolder,foldername)
% filetosave=[pfolder,foldername,'/'];
% for p = fix(clock)
%     filetosave=[filetosave,num2str(p),'_'];
% end
% savefig([filetosave,'fig1','.fig'])

% N_atom=17;
% epslist=0.02:0.02:0.2;
% thr=1e-4;
% % ED_Xelements_sym_func(N_atom,glist,bc,cplist,model,thr)
% for N_atom=N_atomlist
%     cplist=(1:floor(N_atom/2)).^-3;
%     Floquet_Xelements_sym_func(N_atom,epslist,T,bc,cplist,model,thr)
% end
