%20180413
%time evolve an operator using GCR

N_atom=14;
tol=1e-4;
maxiter=100;
bc='o';% 'p' for periodic boundary condition; 'o' for open boundary condition
g=0.5;
% cplist=[1,1/8,1/27];
cplist=1;
delta=0.05;
N=1; % number of steps
tlist=(1:N)*delta;

Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
o_0=Sigma_z;

H_int=sparse(2^N_atom,2^N_atom);

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

col=1;
H_trans=g*Sigma_z;
H=H_int+H_trans;
% U=sparse(1:2*2^N_atom,1:2*2^N_atom,ones(1,2*2^N_atom))+...
%     [sparse(2^N_atom,2^N_atom),delta*H/2;-delta*H/2,sparse(2^N_atom,2^N_atom)];
% UT=U.';
% U2=sparse(1:2*2^N_atom,1:2*2^N_atom,ones(1,2*2^N_atom))+...
%     [sparse(2^N_atom,2^N_atom),delta*H;-delta*H,sparse(2^N_atom,2^N_atom)];
% U2T=U2.';
% rho=full([o_0;zeros(2^N_atom,2^N_atom)]);
% t1=cputime;
% for t=tlist
%     rho=U*rho;
%     parfor p=1:2^N_atom
%         [rho(:,p),res]=tgcr(UT,rho(:,p),tol,maxiter,U2*rho(:,p));
% %         rho(:,p)=U'\rho(:,p);
%     end
%     rho(1:2^N_atom,:)=rho(1:2^N_atom,:).';
%     rho((1+2^N_atom):end,:)=rho((1+2^N_atom):end,:).';
%     rho=UT*rho;
%     parfor p=1:2^N_atom
%         [rho(:,p),res]=tgcr(U,rho(:,p),tol,maxiter,U2T*rho(:,p));
% %         rho(:,p)=U'\rho(:,p);
%     end
%     rho(1:2^N_atom,:)=rho(1:2^N_atom,:).';
%     rho((1+2^N_atom):end,:)=rho((1+2^N_atom):end,:).';
%     t
% end
% t2=cputime-t1;
% rho=rho(1:2^N_atom,:)+1i*rho((1+2^N_atom):end,:);



t1=cputime;
[V,D]=eig(full(H));
t=tlist(end);
rhoED=V*diag(exp(diag(-1i*D*t)))*V'*o_0*V*diag(exp(diag(1i*D*t)))*V';
t3=cputime-t1;
% rhoED=rhoED./sqrt((sum(sum(abs(rhoED).^2))));

% fprintf('max(max(abs(rho-rhoED))=%f\n',max(max(abs(rho-rhoED))))
% fprintf('nnz(U)/(4^N_atom)=%f\n',nnz(U)/(4^N_atom))

% figure(1)
% hold on
% plot(N_atom,t2,'b*')
% plot(N_atom,t3,'k*')
% figure(2)
% hold on
% plot(delta,max(max(abs(rho-rhoED))),'k*')
% figure(3)
% hold on
% plot(N_atom,nnz(U)/(4^N_atom),'k*')