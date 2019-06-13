%20180109
%compute the ground state using exact diagonalization
TT=cputime;
N_atom=12;
N_rand=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

glist=0.5:0.1:1.5;

H_int=zeros(2^N_atom);

for p=1:N_atom-1
%     H_int=H_int+(-1)*...
%         (Sigma_i{1}{p}*Sigma_i{1}{p+1}-...
%         (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{2}{p}*Sigma_i{2}{p+1})/2);
    H_int=H_int+(-1)*...
        (Sigma_i{1}{p}*Sigma_i{1}{p+1});
%     H_int=H_int+...
%         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
    
end

Mzlist=[];
entropy=[];
Mylist=[];
Mxlist=[];

for g=glist
    H_trans=zeros(2^N_atom);
    for p=1:N_atom
        H_trans=H_trans+g*Sigma_i{3}{p};
    end

    H=H_trans+H_int;
    
    [V,D]=eig(H);

    [~,temp]=min(diag(D));
    gs=V(:,temp);
    Mzlist=[Mzlist,gs'*Sigma_z*gs];
    Gs=reshape(gs,[2^(N_atom/2),2^(N_atom/2)]);
    s=svd(Gs);
    entropy=[entropy,-2*log(s')*(s.^2)]; % This may not be correct!!!!
%     Mylist=[Mylist,gs'*Sigma_y*gs];
%     Mxlist=[Mxlist,gs'*Sigma_x*gs];
    g
end

cputime-TT

figure(1)
subplot(4,1,1)
hold on
plot(glist,Mzlist)

subplot(4,1,2)
hold on
plot(glist,entropy)

subplot(4,1,3)
hold on
plot(glist(1:length(glist)-1),Mzlist(2:length(Mzlist))-Mzlist(1:length(Mzlist)-1))

subplot(4,1,4)
hold on
plot(glist(1:length(glist)-1),entropy(2:length(entropy))-entropy(1:length(entropy)-1))