%20180114
%Eigenvalue statics for disorder dipolar model
N_atom=9;
N_rand=100;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,1]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

wlist=logspace(-1,2,20);

rho_0=Sigma_x;
V0=Sigma_z;

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
r=zeros(length(wlist),N_rand);

for p=1:N_atom-1
    H_int=H_int+(-1)*...
        (Sigma_i{1}{p}*Sigma_i{1}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{2}{p}*Sigma_i{2}{p+1})/2);
    %     H_int=H_int+(-1)*...
    %         (Sigma_i{2}{p}*Sigma_i{2}{p+1});
    %     H_int=H_int+...
    %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
end
% H_int=H_int+(-1)*...
%         (Sigma_i{1}{N_atom}*Sigma_i{1}{1}-...
%         (Sigma_i{3}{N_atom}*Sigma_i{3}{1}+Sigma_i{2}{N_atom}*Sigma_i{2}{1})/2);


for pp=1:N_rand
    Disorder=rand(N_atom,1)*2-1;
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(p)*(Sigma_i{1}{p}+Sigma_i{2}{p}+Sigma_i{3}{p})/sqrt(3);
    end
    col=1;
    for w=wlist
        row=1;
        H=H_int+w*H_rand;
        D=eig(H);
        D=sort(D);
        delta=D(2:2^N_atom)-D(1:2^N_atom-1);
        
        for j=1:2^N_atom-2
            r(col,pp)=r(col,pp)+min([delta(j),delta(j+1)])/max([delta(j),delta(j+1)])/(2^N_atom-2);
        end
        col=col+1;
    end
    
end
r_ave=real(sum(r,2))'/N_rand;

if N_rand==1
    figure(2)
    plot(wlist,r_ave/(2*N_atom*2^N_atom))
else
    stdev=std(real(r),1,2)'/N_rand^0.5;
    figure(1)
    hold on
    errorbar(wlist.^-1,r_ave,stdev)
end
