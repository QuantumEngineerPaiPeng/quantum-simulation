% 20190418
% calibrate pulse using the fourier amplitude at the center frequency
% can be used when background signal is significant
% created by Pai Peng
if length(TD)>1
    TD=TD(1);
end
td1=size(FID,3);
td2=size(FID,2);
cFID=FID(1:2:end,:,:)+1i*FID(2:2:end,:,:); % construct complex data
spec=fft((cFID),TD/2,1); 
specs=[spec((end/2+1):end,:,:);spec(1:end/2,:,:)];
amp=squeeze(abs(specs(64,:,:)+specs(63,:,:))/2); % amplitude at center frequence
figure
plot(1:td2,amp);
ampp=amp(:,1); % product of column of amp
for p=2:td1
    ampp=ampp.*amp(:,p);
end
figure
plot(1:td2,ampp)
[~,p2]=max(ampp);
p2