% 20180731
% Transfer a operator in fermionic picture to Pauli picture
% Input: matrix of size 2L by 2L
% Output: matrix of size 2^L by 2^L
% Created by Pai Peng
function y=Fermion2Pauli(M)
M=full(M);
L=size(M,1)/2;
Ide=char(ones(1,L).*'I');
y=sparse(2^L,2^L);
for p=1:L % diagonal terms
    temp=Ide;
    temp(p)='e';
    y=y+PauliProd2Mat(temp)*M(p,p);
    
    temp=Ide;
    temp(p)='g';
    y=y+PauliProd2Mat(temp)*M(p+L,p+L);
end
for p=1:2*L
    for q=1:2*L
        if M(p,q)==0
            continue
        end
        p1=ceil(p/L);
        p2=p-(p1-1)*L;
        q1=ceil(q/L);
        q2=q-(q1-1)*L;
        if q2==p2
            continue
        end
        temp=Ide;
        if p1==1
            temp(p2)='+';
        else
            temp(p2)='-';
        end
        if q1==1
            temp(q2)='-';
        else
            temp(q2)='+';
        end
        pre=1;
        if p2<q2
            for pp=(p2+1):(q2-1)
                temp(pp)='z';
            end
            if temp(p2)=='+'
                pre=-1;
            end
        else
            for pp=(q2+1):(p2-1)
                temp(pp)='z';
            end
            if temp(q2)=='-'
                pre=-1;
            end
        end
%         temp
        y=y+PauliProd2Mat(temp)*M(p,q)*pre;
    end
end
