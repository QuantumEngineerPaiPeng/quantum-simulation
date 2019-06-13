% 20180722
% analysis the Tr(Z(t)Z) data, which shows significant oscillation, so this
% program do the average and then fit.
plist=[4,4,5,6,8,8,10,11,14,20,30];
% plist=ones(1,size(B,2));
% plist=[1,1,5,4,2,2,3,4,5,6,8];
% y=zeros(size(B));
JT=8.18e-3*96*0.2*4;
figure
gamma=zeros(1,size(B,2));
ebar=gamma;
for p=1:size(B,2)
    b1=real(B(:,p));
    b2=zeros(1,length(b1)+1-plist(p));
    for pp=1:(length(b1)+1-plist(p))
        b2(pp)=log(mean(b1(pp:(pp+plist(p)-1))));
    end
%     y(:,p)=b2;
    x=((p+1)/2+(0:(length(b1)-plist(p)))).';
    x=x*JT;
%     plot(x,b2)
    hold on
    f=fit(x,b2.','poly1');
    plot(f,x,b2)
    coe=coeffvalues(f);
    gamma(p)=coe(1);
    conf=confint(f);
    ebar(p)=(conf(2,1)-conf(1,1))/4;
end
figure
hold on
errorbar(XX(1:size(B,2)),-fliplr(gamma),fliplr(ebar))

