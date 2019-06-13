% 20180505
% s=J*s1+J^2*s2/2+J^3*s3/6+...
% find the terms in s^n that is on the order of J^m
function y=findSH(n,m,smax)
y=[];
if n==m
    y=ones(1,n);
else
    if m<n
        y=[];
    else
        if n==1
            if m>smax
                y=[];
            else
                y=m;
            end
        else
            for p=1:min(m,smax)
                temp=findSH(n-1,m-p,smax);
                if ~isempty(temp)
                    y=[y;[p*ones(size(temp,1),1),temp]];
                end
            end
        end
    end
end