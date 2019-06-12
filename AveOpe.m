% 20180711
% compute average operator for a given time points
% arguments: rho, H and tlist
% created by Pai Peng
function y=AveOpe(rho, H, tlist)
if ~IsSameSym({rho,H})
    error('The system size or symmetry class of the input operators are not the same.')
end
if isempty(H.eigsys)
    H.diagonalize
end

rhoV=cell(1,length(rho.matrix));
rhoVbar=cell(1,length(rho.matrix));
rhobar=cell(1,length(rho.matrix));

for p=1:length(rho.matrix)
    rhoV{p}=H.eigsys{p}.V'*rho.matrix{p}*H.eigsys{p}.V;
    rhoVbar{p}=zeros(size(rho.matrix{p}));
end

% y=zeros(length(tlist),length(rho.matrix));
for q=1:length(tlist)
    t=tlist(q);
    for p=1:length(rho.matrix)
        Ut=sparse(diag(exp(-1i*H.eigsys{p}.D*t)));
        rhoVt=Ut*rhoV{p}*Ut';
        rhoVbar{p}=rhoVbar{p}+rhoVt/length(tlist);
    end
end

for p=1:length(rho.matrix)
    rhobar{p}=H.eigsys{p}.V*rhoVbar{p}*H.eigsys{p}.V';
end

y=copy(rho);
y.matrix=rhobar;
y.eigsys=[];