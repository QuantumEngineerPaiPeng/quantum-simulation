% 20180701
% compute two point correlation Tr(rho(t)*O) for an echo sequence with
% forward evolution Hamiltonian H1 and backward evolution Hamiltonian H2
% arguments: rho, H1, H2, O and tlist
% created by Pai Peng
function y=echoTPC(rho, H1, H2, O, tlist)
if ~IsSameSym({rho,H1,H2,O})
    error('The system size or symmetry class of the input operators are not the same.')
end
if isempty(H1.eigsys)
    H1.diagonalize
end

if isempty(H2.eigsys)
    H2.diagonalize
end

rhoV=cell(1,length(rho.matrix));
OV=cell(1,length(rho.matrix));
TM=cell(1,length(rho.matrix));

for p=1:length(rho.matrix)
    rhoV{p}=H1.eigsys{p}.V'*rho.matrix{p}*H1.eigsys{p}.V; % V1'*rho*V1
    OV{p}=H2.eigsys{p}.V'*O.matrix{p}*H2.eigsys{p}.V; % V2'*O*V2
    TM{p}=H1.eigsys{p}.V'*H2.eigsys{p}.V; % V1'*V2
end

y=zeros(length(tlist),length(rho.matrix));
for q=1:length(tlist)
    t=tlist(q);
    for p=1:length(rho.matrix)
        U1t=sparse(diag(exp(-1i*H1.eigsys{p}.D*t)));
        U2t=sparse(diag(exp(-1i*H2.eigsys{p}.D*t)));
%         rhoVt=U1t*rhoV{p}*U1t';
%         OVt=U2t'*OV{p}*U2t;
%         y(q,p)=sum(sum((TM{p}*rhoVt*TM{p}').*(OVt.')));
        y(q,p)=trace(U2t*TM{p}'*U1t*rhoV{p}*U1t'*TM{p}*U2t'*OV{p});
    end
end

y=sum(y,2);