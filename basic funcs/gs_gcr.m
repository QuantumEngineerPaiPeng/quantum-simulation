% 20181027
% compute the ground state using gcr
% Written using operator class

function gs=gs_gcr(H)
% DisorderDir=[0,0,1]; %direction of the disorder field
% DisorderDir=DisorderDir/norm(DisorderDir);
% 
% N_atom=10;
% H_int=OperatorClass(N_atom,'dipolar',-1,'p',(1:N_atom-1).^-3,1);
% rng(1)
% randmat=rand(1,N_atom);
% Disorder=DisorderDir.'*(randmat*2-1);
% H_rand=randOpe(N_atom,Disorder);
% H=H_int+1*H_rand;

if isa(H,'OperatorClass')
    if ~isempty(H.sym)
        error('Operator to calculate ground state must not be symmetrized.')
    end
    M=H.matrix{1};
else
    M=H;
end

tol=1e-8;
maxiter=100;


phi0=zeros(length(M),1);
phi0(length(M))=1;


ns=1;
U=0.01*(M)+sparse(1:length(M),1:length(M),ones(1,length(M)));

phi=[phi0,tgcr(U,phi0,tol,maxiter,phi0)];
phi(:,2)=phi(:,2)/norm(phi(:,2));

while norm(phi(:,1)-phi(:,2))>1e-6
    phi(:,1)=phi(:,2);
    phi(:,2)=tgcr(U,phi(:,2),tol,maxiter,phi(:,2));
    phi(:,2)=phi(:,2)/norm(phi(:,2));
    ns=ns+1;
end
gs=phi(:,2);

