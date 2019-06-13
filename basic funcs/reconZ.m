% 20180920
% Inverse process of DecomZ
% Input: 2^L array
% Created by Pai Peng
function y= reconZ(x)
sl=1; % sector length
L=round(log2(length(x)));
for p=1:L
    x1=reshape(x,sl,[]);
    x2=zeros(size(x1));
    for q=1:size(x1,2)/2
        x2(:,2*q-1)=x1(:,2*q-1)+x1(:,2*q);
        x2(:,2*q)=x1(:,2*q-1)-x1(:,2*q);
    end
    x=reshape(x2,1,[]);
    sl=sl*2;
end
y=x.';