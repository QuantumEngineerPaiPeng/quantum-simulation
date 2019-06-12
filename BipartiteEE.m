% 20180812
% Get the bipartite entanglement entropy
% Input: 
% rho, given in OperatorClass
% n, partition at n, n+1
% Created by Pai Peng
function y=BipartiteEE(rho,n)
% rho.normalize;
if abs(rho.trace-1)>1e-12
    rho.matrix{1}=rho.matrix{1}/rho.trace;
end
if ~isempty(rho.sym)
    error('Operator must not be symmetrized')
end
% matr=rho.matrix{1};
% for p=1:n
%     matr=matr(1:end/2,1:end/2)+matr((end/2+1):end,(end/2+1):end);
% end
matr=zeros(2^(rho.L-n));
for p=1:2^n
    matr=matr+rho.matrix{1}(((p-1)*2^(rho.L-n)+1):p*2^(rho.L-n),((p-1)*2^(rho.L-n)+1):p*2^(rho.L-n));
end
D=eig(full(matr));
dd=[];
for p=1:length(D)
    if abs(D(p))>1e-13
        dd=[dd,D(p)];
    end
end
if min(dd)<0
    error('rho is not positive definite')
end
y=-dd*log(dd.');