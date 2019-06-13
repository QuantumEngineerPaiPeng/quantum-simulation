% OTOC_t_g_2gs_sym_func(18,0.4:0.2:3,1:1:15,'p',[1,1/8,1/27],'dipolar')

% Natomlist=8;
% H0='xx';
% bc='o';
% 
% tic
% for N_atom=Natomlist
%     cplist=(1:N_atom-1).^(-3);
%     prethermalH_func(N_atom,N_atom,bc,cplist,H0)
% end
% toc

% N_atom=11;
% order=11;
% epslist=0.25:0.25:3;
% pre_r_order_func(N_atom,order,epslist)

% tic
% prexx_r_order_func(8,8,0.25:0.25:1)
% toc

% N_atom=10;
% order=N_atom;
% J=0.25;
% glist=logspace(-0.8,0,5);
% bc='o';
% cplist=(1:N_atom-1).^(-3);
% for g=glist
%     prethermalHxx_func(N_atom,g,J,order,bc,cplist)
% end

% f1_g

% OTOC_t_g_sym_any

% s_locality

t_disorder