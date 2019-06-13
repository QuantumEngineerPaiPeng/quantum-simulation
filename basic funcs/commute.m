% 20180327
% for given input a,b='x','y', or 'z', representing Pauli matrices return [a,b]
function [coe,op]=commute(a,b)
if or(a=='1',b=='1')
    coe=0;
    op='';
else
    a=a-'x';
    b=b-'x';
    if a==b
        coe=0;
        op='';
    else
        op=char('x'+3-a-b);
        if mod(a+1,3)==b
            coe=2*1i;
        else
            coe=-2*1i;
        end
    end
end
end