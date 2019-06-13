% 20180830
% Decompose the density matrix according to both Hamming distance and
% correlation length, and get the weight in each sector.
% input
% O: the matrix or OperatorClass
% opt1: 'sum' or 'mean', specifies taking the sum or mean in each sector
% opt2: 'square' or 'abs', specifies taking the squared or absolute value
% of each element
% Created by Pai Peng
function z=HammingCLWeight(O,opt1,opt2)
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

z=zeros(L+1,L+1);
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/CLcell%d_%d.mat',L,d);
if exist(filename,'file')
    load(filename)
else
    y=CLcell(L,d);
end
CL=y;

filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/Hammingcell%d_%d.mat',L,d);
if exist(filename,'file')
    load(filename)
else
    y=Hammingcell(L,d);
end
Hamming=y;

for p=1:L+1 % index for Hamming distance
    for q=p:L+1 % index for correlation length
        if strcmp(opt2,'abs')
            temp=abs(amp(intersect(Hamming{p},CL{q})));
            temp(temp<1e-13)=[];
        else
            if strcmp(opt2,'square')
                temp=abs(amp(intersect(Hamming{p},CL{q})));
                temp(temp<1e-13)=[];
                temp=temp.^2;
            else
                error('3rd argument must be ''abs'' or ''square''.')
            end
        end
        if isempty(temp)
            z(p,q)=0;
        else
            z(p,q)=eval([opt1,'(temp)']);
        end
    end
end