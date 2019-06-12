% 20180327
% for a given operator o (sum of product operators such as 1*'111'+2*'xyz'),
% and a Hamiltonian H (sum of product operators such as 1*'xx'+0.5*'z')
% return the [o,H] in a cell and the corresponding coefficient in a
% list

function [coe,result]=commuteOH(coe0,o,coeH,H)
% o={'111xy111','1111yx11'};
% coe0=[1,1];
%
% H={'xx','z'};
% coeH=[1,1];

coe=[];
result={};
for p=1:length(o)
    for pp=1:length(H)
        [c1,o1]=commute_prod(o{p},H{pp});
        for ppp=1:length(o1)
            result{end+1}=o1{ppp};
            coe(end+1)=coe0(p)*coeH(pp)*c1(ppp);
        end
    end
end

% combine duplicate operators
for p=1:length(result)
    for pp=p+1:length(result)
        if strcmp(result{p},result{pp})
            coe(p)=coe(p)+coe(pp);
%             result(pp)=[];
            coe(pp)=0;
        end
    end
end

% delete operators with 0 coefficient
p=1;
while p<=length(result)
    if abs(coe(p))<1e-10
        result(p)=[];
        coe(p)=[];
    else
        p=p+1;
    end
end