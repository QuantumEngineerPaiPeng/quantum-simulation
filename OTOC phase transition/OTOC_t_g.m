%20180108
%compute OTOC v.s g for different t
tic
N_atom=10;
N_rand=1;
bc='p';% 'p' for periodic boundary condition; 'o' for open boundary condition
cplist=[1,1/8,1/27];
% cplist=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

%tlist=(1.88:1.88:8)./4;
tlist=1:0.5:3;
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];
glist=.2:0.2:3;

rho_0=Sigma_z;
V0=Sigma_x;

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
C=zeros(length(tlist),length(glist),N_rand);
for pp=1:length(cplist)
    for p=1:N_atom
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else 
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
%         H_int=H_int+(-1)*cplist(pp)*...
%             (Sigma_i{1}{p}*Sigma_i{1}{q});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end
% H_int=H_int+(-1)*...
%     (Sigma_i{1}{1}*Sigma_i{1}{N_atom}-...
%     (Sigma_i{3}{1}*Sigma_i{3}{N_atom}+Sigma_i{2}{1}*Sigma_i{2}{N_atom})/2);
%hh=waitbar(0);
for pp=1:N_rand
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(rand-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(rand-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(rand-0.5)*Sigma_i{3}{p};
    end
    col=1;
    for g=glist
        H_trans=g*Sigma_z;
        H=H_int+H_trans+H_rand;
        [V,D]=eig(H);
        row=1;
        for t=tlist
            
            rhot=V*diag(exp(diag(1i*D*t)))*V'*rho_0*V*diag(exp(diag(-1i*D*t)))*V';
            C(row,col,pp)=-trace((rhot*V0-V0*rhot)^2);
            row=row+1;
        end
        %hh=waitbar(col/length(glist));
        col=col+1;

        g
    end
    
end

C_ave=real(sum(C,3))'/N_rand/(2*N_atom*2^N_atom);
[maxC,temp]=max(C_ave,[],1);
maxg=glist(temp);
% clearvars('H_int','H_rand','H','H_trans')

figure(1)
hold on
plot(tlist,maxg)
figure(2)
hold on
plot(tlist,maxC)
if N_rand==1
    figure(3)
    hold on
    plot(glist,(C_ave)/2,'r','LineWidth',2)
else
    stdev=std(real(C),1,3)'/N_rand^0.5;
    figure(3)
    hold on
    for p=1:N_t
        errorbar(tlist,C_ave(:,p),stdev(:,p))
    end
end
toc