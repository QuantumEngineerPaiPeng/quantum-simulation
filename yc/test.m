%%
tau_list_all = tau_list;
r_list_all = r_list;
dr_list_all = dr_list;

%%
%tau_list_all = [tau_list_all, tau_list];
r_list_all = cat(2, r_list(:,:,1:54), r_list_all);
dr_list_all = cat(2, dr_list(:,:,1:54), dr_list_all);

%% 
tau_list = tau_list_all;
r_list = r_list_all;
dr_list = dr_list_all;

%%
[u, v, w] = [1, 2, 3]
