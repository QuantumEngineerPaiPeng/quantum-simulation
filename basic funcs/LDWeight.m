% 20180827
% Calculate the localization distance weight of an OperatorClass or matrix
% e.g. IIXYXZII, IXIIIXII are of distance 3 to 3rd site
% input:
% O: operator to evalute
% s: the site to which the distance is calculated
% opt1: 'sum' or 'mean', specifies taking the sum or mean in each sector
% opt2: 'square' or 'abs', specifies taking the squared or absolute value
% Created by Pai Peng
function z=LDWeight(O,s,opt1,opt2)
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
if isvector(matr)
    d=2;
    amp=DecomZ(matr);
else
    d=4;
    amp=DecomPauli(matr);
end
y=LDcell(L,d,s);
z=zeros(L,1);
for p=1:L
    if strcmp(opt2,'abs')
        temp=abs(amp(y{p}));
        temp(temp<1e-13)=[];
    else
        if strcmp(opt2,'square')
            temp=abs(amp(y{p}));
            temp(temp<1e-13)=[];
            temp=temp.^2;
        else
            error('3rd argument must be ''abs'' or ''square''.')
        end
    end
    if isempty(temp)
        z(p)=0;
    else
        z(p)=eval([opt1,'(temp)']);
    end
end