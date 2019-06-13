% 20180711
% Analysis for averaged OTOC, pp e.g. xw_dqxy_mblotoavg1
tic
N_phi=8;
N_t1=6;
N_t2=6;
otoc=mblOTOCavg_func(B,N_phi,N_t1,N_t2);
% N_phi=6;
% N_t=12;
% N_cycle=4;
% 
% otoc=OTOC_exp_func(B,N_phi,N_cycle,N_t);
VAR=zeros(size(otoc));
for p=1:size(B,1)
    for q=1:size(B,2)
        errB=zeros(size(B));
        errB(p,q)=err(p,q);
        VAR=VAR+(mblOTOCavg_func(B+errB,N_phi,N_t1,N_t2)-otoc).^2;
%         VAR=VAR+(OTOC_exp_func(B+errB,N_phi,N_cycle,N_t)-otoc).^2;
    end
end
ebar=sqrt(VAR);
toc