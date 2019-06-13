% 20180808
% fits the Tr(Y(t)Y) data (1613:1623) to damped oscllation.
% data stored in B
% Created by Pai Peng

figure
gamma=zeros(1,size(B,2));
r2=gamma;
JT=8.18e-3*96*0.2*4;
x=((1:64)*JT).';
ebar=gamma;
for p=1:size(B,2)
    b1=real(B(:,p));

    fun=@(a,b,c,gamma,x) a*exp(-gamma*x).*cos(b*x+c);

    hold on
%     if p==1
        [fitted,gof]=fit(x,b1,fun,'StartPoint',[8e3,0,0.01,0.1],'Lower',[6e3,0,0,0],'Upper',[1e4,10,2*pi,Inf]);
%     else
%         [fitted,gof]=fit(x,b1,fun,'StartPoint',coeffvalues(fitted));
%     end
    plot(fitted,x,b1)
    coe=coeffvalues(fitted);
    gamma(p)=coe(4);
    conf=confint(fitted);
    ebar(p)=(conf(2,4)-conf(1,4))/4;
    r2(p)=gof.rsquare;

end
figure
hold on
errorbar(XX(1:size(B,2)).^-1,fliplr(gamma),fliplr(ebar))

