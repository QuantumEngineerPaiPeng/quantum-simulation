% 20180620
% compute two point correlation Tr(rho(t)*O)
% arguments: rho, H, O and tlist
% created by Pai Peng
function y=TPC(rho, H, O, tlist)
if ~IsSameSym({rho,H,O})
    error('The system size or symmetry class of the input operators are not the same.')
end
if isempty(H.eigsys)
    H.diagonalize
end

rhoV=cell(1,length(rho.matrix));
OV=cell(1,length(rho.matrix));

for p=1:length(rho.matrix)
    rhoV{p}=H.eigsys{p}.V'*rho.matrix{p}*H.eigsys{p}.V;
    OV{p}=H.eigsys{p}.V'*O.matrix{p}*H.eigsys{p}.V;
end

y=zeros(length(tlist),length(rho.matrix));
for q=1:length(tlist)
    t=tlist(q);
    for p=1:length(rho.matrix)
        Ut=sparse(diag(exp(-1i*H.eigsys{p}.D*t)));
        rhoVt=Ut*rhoV{p}*Ut';
        y(q,p)=sum(sum(rhoVt.*(OV{p}.')));
    end
end

y=sum(y,2);