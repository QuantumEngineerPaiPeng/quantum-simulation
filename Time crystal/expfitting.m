% 20180604
% exponential fit to DTC data
Bb=squeeze(abs(B));
frame=(mod(1:size(Bb,1),2)-0.5)*2;
% Bb=Bb.*frame';
iu10=4;
x=0:(size(Bb,1)-1);
x=1+iu10*x;
varname='perturbation';
% perf_pi=1.95;
% Xlist=(2.15:-0.05:1.7)-perf_pi;
% perf_pi=2;
% Xlist=(1.5:0.05:2.45)-perf_pi;
% varname='L';
% Xlist=10:17;
% Xlist=epslist(1:10);
% varname='period';
% Xlist=[36,50,60,70,80,100];
figure
hold on
r2list=zeros(1,size(Bb,2));
paralist=zeros(2,size(Bb,2));
errlist=zeros(2,size(Bb,2));
fitlist=cell(1,size(Bb,2));
goflist=cell(1,size(Bb,2));
for p=1:size(Bb,2)
    if p>1
        [fitted,gof] = fit(x.',Bb(:,p),'exp1','StartPoint',paralist(:,p-1));
    else
        [fitted,gof] = fit(x.',Bb(:,p),'exp1');
    end
    fitlist{p}=fitted;
    goflist{p}=gof;
    r2list(p)=gof.rsquare;
    paralist(:,p)=coeffvalues(fitted).';
    temp=confint(fitted);
    errlist(:,p)=((temp(2,:)-temp(1,:))/2).';
    
    try
        plot(fitted)
        errorbar(x,Bb(:,p),err(:,p))
    catch
        plot(fitted,x,Bb(:,p))
    end
        xlabel('# of Floquet cycles')
        ylabel('|M|')
        ppStyle(30,2,15)

end
figure
plot(Xlist,r2list)
xlabel(varname)
ylabel('R^2')

figure
for p=1:2
    subplot(1,2,p)
    hold on
    try
        errorbar(Xlist,paralist(p,:),errlist(p,:))
    catch
        plot(Xlist,paralist(p,:))
    end
    ppStyle(15,2)
    xlabel(varname)
    ylabel(char('a'-1+p))
    title('fit to a*exp(b*x)')
end