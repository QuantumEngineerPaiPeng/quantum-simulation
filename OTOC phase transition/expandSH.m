% 20180505
% find the nth order expansion of e^s*H*e^-s
% with s=J*s1+J^2*s2/2+J^3*s3/6+...
% see https://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula
function y=expandSH(s,H,order)
if order==0
    y=H;
else
    if isempty(s)
        y=sparse(size(H,1),size(H,2));
    else
        y=sparse(size(H,1),size(H,2));
        ytemp=cell(1,order);
        parfor p=1:order
            list=findSH(p,order,length(s));
            yy=sparse(size(H,1),size(H,2));
            for pp=1:size(list,1)
                temp=H;
                %             p
                %             list(pp,:)
                %             pause
                for ppp=1:p
                    %                 s(:,:,list(pp,ppp))
                    temp=s{list(pp,ppp)}*temp-temp*s{list(pp,ppp)};
                    %                 pause
                end
                yy=yy+temp/factorial(p);
            end
            ytemp{p}=yy;
        end
        for p=1:order
            y=y+ytemp{p};
        end
    end
end