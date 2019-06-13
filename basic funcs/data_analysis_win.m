numStart=572;
numEnd=572;
A=[];
for p=numStart:numEnd
    %filename=sprintf('/Volumes/nmrsu1/expt11/%d/ser',p);
    filename=sprintf('~/Dropbox (MIT)/QEG/NMR_manual/data/%d/ser',p);
    fid=fopen(filename,'r','b');
    A=[A,fread(fid,'int32')];
end
B=reshape(A,256,5);
C=B(1:2:256,:)+1i*B(2:2:256,:);
% frame=mod((1:128),2)*2-1;
% frameCr=frame'.*Cr;
figure(1)
hold on
plot(1:128,abs(C))