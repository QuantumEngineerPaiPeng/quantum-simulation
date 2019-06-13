%20180210
%load and analysis data from .m files
epslist=[];
gammalist=[];
errlist=[];
cf=[];
filelist=dir(fullfile('DTC_0.5*'));
for p=1:length(filelist)
    file=filelist(p);
    load(file.name)
    %if ~exist('tlist','var')
        tlist=(0:length(MExp)-1)*0.5;
    %end
    frame=1:length(MExp);
    frame=mod(frame,2)*2-1;
    fitp=fit(tlist.',(frame.*MExp).','exp1')
    filename=file.name;
    temp=strread(filename(1:length(filename)-4),'%s','delimiter','_');
    if ~exist('epsilon','var')
        epsilon=str2num(temp{3});
    end
    epslist=[epslist,epsilon];
    coeff=coeffvalues(fitp);
    gammalist=[gammalist,-coeff(2)];
    confid=confint(fitp);
    errlist=[errlist,[confid(2,2)-coeff(2);coeff(2)-confid(1,2)]];
    L=2^floor(log2(length(MExp)));
    famp=fft(MExp(length(MExp)-L+1:length(MExp)),L);
    %famp=fft(MExp,L);
    cf=[cf,abs(famp(L/2+1))^2/(sum(abs(famp).^2))];
    clearvars('tlist')
    clearvars('epsilon')
end
figure(1)
%errorbar(epslist,gammalist,errlist(1,:),errlist(2,:))
plot(epslist,gammalist)
figure(2)
plot(epslist,cf)