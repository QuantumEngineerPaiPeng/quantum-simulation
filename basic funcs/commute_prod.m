% 20180327
% for a given product operator o, such as '111xyz1yz111', and a product Hamiltonian H, such as
% 'xx', return the [o,H] in a cell and the corresponding coefficient in a
% list
function [coe,result]=commute_prod(o,H)

lH=length(H);
lo=length(o);
result={};
coe=[];
for p=0:(lo-lH)
    for pp=1:lH
        [c1,o1]=commute(o(p+pp),H(pp));
        if c1~=0
            temp=o;
            temp(p+pp)=o1;
            for ppp=1:lH
                if ppp<pp
                    [c2,o2]=timep(o(p+ppp),H(ppp));
                    temp(p+ppp)=o2;
                    c1=c1*c2;
                else
                    if ppp>pp
                        [c2,o2]=timep(H(ppp),o(p+ppp));
                        temp(p+ppp)=o2;
                        c1=c1*c2;
                    end
                end
            end
            coe(end+1)=c1;
            result{end+1}=temp;
        end
    end
end
end