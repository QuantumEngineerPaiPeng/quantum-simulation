%20180402
%compute OTOC with local operators v.s position for different t
tic
N_atom=9;
N_rand=1;
bc='o';% 'p' for periodic boundary condition; 'o' for open boundary condition
o1=[5,3]; % time evolved operator, 1st para: position, 2nd para: Pauli
o2=1; % constant operator, with only one para being the Pauli
g=2;
% cplist=[1,1/8,1/27];
cplist=1;

Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

%tlist=(1.88:1.88:8)./4;
tlist=linspace(0,10,30);
[tt,pos]=meshgrid(tlist,1:N_atom);


rho_0=Sigma_i{o1(2)}{o1(1)};

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
C=zeros(N_atom,length(tlist),N_rand);
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
%                 H_int=H_int+(-1)*cplist(pp)*...
%                     (Sigma_i{1}{p}*Sigma_i{1}{q}-...
%                     (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q});
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
    H_trans=g*Sigma_z;
    H=H_int+H_trans+H_rand;
    [V,D]=eig(H);
    
    for t=tlist
        row=1;
        rhot=V*diag(exp(diag(1i*D*t)))*V'*rho_0*V*diag(exp(diag(-1i*D*t)))*V';
        for p=1:N_atom
            V0=Sigma_i{o2}{p};
            C(row,col,pp)=-trace((rhot*V0-V0*rhot)^2);
            row=row+1;
        end
        col=col+1;
    end
end

C_ave=real(sum(C,3))/N_rand/(4*2^N_atom);
figure
h=pcolor(pos,tt,C_ave);
set(h, 'EdgeColor', 'none');
colorbar

% shading interp
toc