% 20180808
% Get the time whe the Tr(Y(t)Y) data (1613:1623) decays to thr*initial
% value
% data stored in B
% Created by Pai Peng

thr=1/e;
tt=zeros(1,size(B,2));
r2=gamma;
JT=8.18e-3*96*0.2*4;
x=((1:64)*JT).';
ebar=gamma;
for p=1:size(B,2)
    b1=real(B(:,p));
    for q=2:size(B,1)
        if b1(q)<b1(1)*thr
            tt(p)=q+(b1(q)-b1(1)*thr)/(-(b1(q)-b1(1)*thr)+b1(q-1)-b1(1)*thr);
            q
            break
        end
    end
end
tt=tt*JT;
figure
hold on
plot(XX(1:size(B,2)).^-1,fliplr(tt.^-1))

