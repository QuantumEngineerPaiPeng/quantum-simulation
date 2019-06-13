% 20180911
% Calculate the localization distance weight of an matrix in Fermion
% picture
% input:
% O: operator to evalute
% s: the site to which the distance is calculated
% opt1: 'sum' or 'mean', specifies taking the sum or mean in each sector
% opt2: 'square' or 'abs', specifies taking the squared or absolute value
% Created by Pai Peng
% See also LDWeight
function z=LDWeightFermion(O,s,opt1,opt2)
L=length(O)/2;
O1=O(1:L,1:L);
O2=O(1:L,(L+1):end);
z=zeros(L,1);
count=zeros(L,1);
for p=1:L
    for q=p:L
        if abs(O1(p,q))>1e-13
            d=max(abs(p-s),abs(q-s));
            if strcmp(opt2,'square')
                temp=abs(O1(p,q))^2;
            else
                temp=abs(O1(p,q));
            end
            if p~=q
                temp=temp*4;
                count(d+1)=count(d+1)+4;
            else
                temp=temp*2;
                count(d+1)=count(d+1)+2;
            end
            z(d+1)=z(d+1)+temp;
        end
        if abs(O2(p,q))>1e-13
            d=max(abs(p-s),abs(q-s));
            if strcmp(opt2,'square')
                temp=abs(O2(p,q))^2;
            else
                temp=abs(O2(p,q));
            end
            if p~=q
                temp=temp*4;
                count(d+1)=count(d+1)+4;
            else
                temp=temp*2;
                count(d+1)=count(d+1)+2;
            end
            z(d+1)=z(d+1)+temp;
        end
    end
end

if strcmp(opt1,'mean')
    for p=1:L
        if count(p)>0
            z(p)=z(p)/count(p);
        else
            z(p)=0;
        end
    end
end