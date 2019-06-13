%20180402
%time evolve an operator using forward Euler

N_atom=9;

bc='o';% 'p' for periodic boundary condition; 'o' for open boundary condition
g=0.5;
% cplist=[1,1/8,1/27];
cplist=1;
delta=1e-4;
beta=1e-4;
N=100; % number of steps
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
o0=Sigma_z;
rho_0=eye(2^N_atom)-o0;
% rho_0h=sparse(sqrt(rho_0));
rho_0h=eye(2^N_atom)-beta*o0/2;

H_int=zeros(2^N_atom);

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
U=sparse(eye(2^N_atom)-1i*delta*H-1/2*delta^2*H^2);
% U=sparse(eye(2^N_atom)-1i*delta*H);
t1=cputime;
rhoh=rho_0h;
for t=tlist
    rhoh=U*rhoh;
end
t2=cputime-t1;

rho=(rhoh-sparse(eye(2^N_atom)))/(-beta/2);
% rho=-rho./sqrt((sum(sum(abs(rho).^2))));
t1=cputime;
[V,D]=eig(H);
rhoED=V*diag(exp(diag(-1i*D*t)))*V'*o0*V*diag(exp(diag(1i*D*t)))*V';
t3=cputime-t1;
% rhoED=rhoED./sqrt((sum(sum(abs(rhoED).^2))));

% fprintf('max(max(abs(rho-rhoED))=%f\n',max(max(abs(rho-rhoED))))
% fprintf('nnz(U)/(4^N_atom)=%f\n',nnz(U)/(4^N_atom))

figure(1)
hold on
plot(N_atom,t2,'r*')
plot(N_atom,t3,'k*')
figure(2)
hold on
plot(N_atom,max(max(abs(rho-rhoED))),'k*')
figure(3)
hold on
plot(N_atom,nnz(U)/(4^N_atom),'k*')