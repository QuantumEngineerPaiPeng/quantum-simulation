% 20180421
% calibrate frequency
FF=reshape(FID,TD,[]);
FF=FF(1:2:end,:)+1i*FF(2:2:end,:);
fourier=fft(FF,size(FF,1),1);
fourier=[fourier(round(size(fourier,1)/2+1):end,:);fourier(1:round(size(fourier,1)/2),:)];
figure
for p=1:size(fourier,2)
    subplot(3,4,p)
    hold on
    plot(1:size(fourier,1),abs(fourier(:,p)))
    plot(size(fourier,1):-1:1,abs(fourier(:,p)))
end