N_atomlist=10:2:12;
Nlist=1:100:1e3;
epslist=0:0.5:1;
T=0.5;
bc='p';
cplist=[1,1/8,1/27];
model='dipolar';
figure(1)
hold on
for N_atom=N_atomlist
OTOC_t_g_2gs_sym_func(N_atom,Nlist,epslist,T,bc,cplist,model)

end
pfolder='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
foldername=mfilename;
mkdir(pfolder,foldername)
filetosave=[pfolder,foldername,'/'];
for p = fix(clock)
    filetosave=[filetosave,num2str(p),'_'];
end
savefig([filetosave,'fig1','.fig'])