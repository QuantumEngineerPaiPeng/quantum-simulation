% 20180817
% Calculate the corration length weight of an OperatorClass or matrix
% e.g. IIXYXZII, IXIIXIII are of length 4
% Created by Pai Peng
% Replaced by CLWeight on 20180829
function y=CLWeight_Old(O)
% H_int=OperatorClass(8,'Ising',-1,'o',(1:3).^-3,1);
if isa(O,'OperatorClass')
    if ~isempty(O.sym)
        error('Operator must not be symmetrized')
    end
    L=O.L;
    matr=O.matrix{1};
else
    matr=O;
    L=round(log2(length(matr)));
end
amp=DecomPauli(matr);
y=zeros(L+1,1);
for p=0:(length(amp)-1)
    if abs(amp(p+1))>1e-14
        ps=dec2base(p,4);
        while ps(end)=='0'
            ps(end)=[];
            if isempty(ps)
                break;
            end
        end
        
        y(length(ps)+1)=y(length(ps)+1)+amp(p+1)*amp(p+1)';
    end
end