%20180402
%compute the butterfly velocity OTOC with local operators v.s g
tic
N_atom=8;
N_rand=1;
bc='o';% 'p' for periodic boundary condition; 'o' for open boundary condition
o1=[2,1]; % time evolved operator, 1st para being the Pauli, 2nd being position
o2=[3,N_atom]; % constant operator, 1st para being the Pauli, 2nd being position
glist=linspace(0,3,30);
thr=0.4;
dg=1e-3;
gdg=reshape([glist;glist+dg],1,2*length(glist));
tlist=zeros(size(glist));
% cplist=[1,1/8,1/27];
cplist=1;

Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,0,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]

rho_0=Sigma_i{o1(1)}{o1(2)};
V0=Sigma_i{o2(1)}{o2(2)};

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
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
%                 H_int=H_int+(-1)*cplist(pp)*...
%                     (Sigma_i{1}{p}*Sigma_i{1}{q});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end

for pp=1:N_rand
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(rand-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(rand-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(rand-0.5)*Sigma_i{3}{p};
    end
    col=1;
    for g=gdg
        H_trans=g*Sigma_z;
        H=H_int+H_trans+H_rand;
        [V,D]=eig(H);
        t0=0;
        t1=10;
        while t1-t0>1e-5
            t=(t0+t1)/2;
            rhot=V*diag(exp(diag(1i*D*t)))*V'*rho_0*V*diag(exp(diag(-1i*D*t)))*V';
            C=-real(trace((rhot*V0-V0*rhot)^2))/(4*2^N_atom);
            if C<thr
                t0=t;
            else
                t1=t;
            end
        end
        tlist(col)=(t0+t1)/2;
        col=col+1;
    end
end
vB=(o2(2)-o1(2))*tlist.^(-1);
figure(1)
hold on
plot(glist,vB(1:2:end),'DisplayName',[num2str(o1(1)),',',num2str(o2(1))])
xlabel('g')
ylabel('v_B')
ppStyle(30,2)
legend()

figure(2)
hold on
plot(glist+dg/2,(vB(2:2:end)-vB(1:2:end))/dg,'DisplayName',[num2str(o1(1)),',',num2str(o2(1))])
xlabel('g')
ylabel('dv_B/dg')
legend()
ppStyle(30,2)
% shading interp
toc