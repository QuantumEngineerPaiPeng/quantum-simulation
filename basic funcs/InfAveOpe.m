% 20180803
% compute infiite average operator by keeping only the diagonal elements in
% the eigebasis
% arguments: rho, H
% created by Pai Peng
function y=InfAveOpe(rho, H)
if ~IsSameSym({rho,H})
    error('The system size or symmetry class of the input operators are not the same.')
end

% if ~isempty(H.sym)
%     error('Operators must not be symmetrized')
% end

if isempty(H.eigsys)
    H.diagonalize()
end

y=copy(H);
for q=1:length(H.matrix)
    temp=zeros(size(H.matrix{q},1),1);
    for qq=1:size(H.eigsys{1}.V,1)
        temp(qq)=H.eigsys{q}.V(:,qq)'*rho.matrix{q}*H.eigsys{q}.V(:,qq);
    end
    y.matrix{q}=H.eigsys{q}.V*sparse(diag(temp))*H.eigsys{q}.V';
    y.eigsys{q}.D=temp;
end



