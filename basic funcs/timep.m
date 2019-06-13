% 20180327
% for given input a,b='x','y', or 'z', representing Pauli matrices return
% ab
function [coe,op]=timep(a,b)
if a=='1'
    coe=1;
    op=b;
else
    if b=='1'
        coe=1;
        op=a;
    else
        a=a-'x';
        b=b-'x';
        if a==b
            coe=1;
            op='1';
        else
            op=char('x'+3-a-b);
            if mod(a+1,3)==b
                coe=1i;
            else
                coe=-1i;
            end
        end
    end
end