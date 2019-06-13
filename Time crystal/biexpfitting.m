% 20180415
% biexponential fit to DTC data
% Bb=squeeze(-real(B));
% Bb=Bb./(Bb(1,:));
Bb=-real(B);
% Bb=real(Bb);
% Bb=squeeze(abs(B));
frame=(mod(1:size(Bb,1),2)-0.5)*2;
% Bb=Bb.*frame';
% Bb=real(Mt(2:end,:)).';
iu10=1;
% xname='# of Floquet cycles';

x=0:(size(Bb,1)-1);
x=1+iu10*x;
x=x*72;
xname='t (us)';
varname='perturbation';
perf_pi=2.06;
Xlist=((2.36:-0.05:1.76)-perf_pi)/perf_pi*pi;

% Xlist=epslist(2:end);
% Xlist=epslist;
% Xlist=-0.4:0.01:-0.06;
% varname='L';
% Xlist=10:17;

% varname='period';
% Xlist=[36,50,60,70,80,100];
figure
hold on
r2list=zeros(1,size(Bb,2));
paralist=zeros(4,size(Bb,2));
errlist=zeros(4,size(Bb,2));
fitlist=cell(1,size(Bb,2));
goflist=cell(1,size(Bb,2));

fun=@(a,b,c,d,x) log(-a*exp(b*x)-c*exp(d*x));

for p=1:size(Bb,2)
%     if p>1
%         [fitted,gof] = fit(x.',Bb(:,p),'exp2','StartPoint',paralist(:,p-1));
%     else
%         [fitted,gof] = fit(x.',Bb(:,p),'exp2');
%     end
    if p>1
        [fitted,gof] = fit(x.',log(-Bb(:,p)),fun,'StartPoint',paralist(:,p-1));
    else
        [fitted,gof] = fit(x.',log(-Bb(:,p)),fun,'StartPoint',[-5000,-0.002,-5000,-5e-4]);
    end
    fitlist{p}=fitted;
    goflist{p}=gof;
    r2list(p)=gof.rsquare;
    paralist(:,p)=coeffvalues(fitted).';
    temp=confint(fitted);
    errlist(:,p)=((temp(2,:)-temp(1,:))/2).';
%     plot(fitted,x,Bb(:,p))
    plot(fitted,x,log(-Bb(:,p)))
%     errorbar(x,Bb(:,p),err(:,p))
    xlabel(xname)
    ylabel('|M|')
    ppStyle(30,2,15)
end
figure(4)
plot(Xlist,r2list)
xlabel(varname)
ylabel('R^2')

figure(6)
for p=1:4
    subplot(2,2,p)
    hold on
    errorbar(Xlist,paralist(p,:),errlist(p,:))
    ppStyle(15,2)
    xlabel(varname)
    ylabel(char('a'-1+p))
    title('fit to a*exp(b*x)+c*exp(d*x)')
end