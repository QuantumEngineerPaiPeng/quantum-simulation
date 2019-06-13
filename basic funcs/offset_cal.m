% 20180618
% offset calibration
if length(TD)>1
    TD=TD(1);
end
normalize=false;
cFID=FID(1:2:end,:)+1i*FID(2:2:end,:); % construct complex data
spec=fft((cFID),TD/2,1); 
specs=[spec((end/2+1):end,:);spec(1:end/2,:)];
specd=abs(spec((end/2+1):end,:))-abs(spec(end/2:-1:1,:));
% specd=
figure
if normalize
    plot(((0:(TD/2-1))-(TD/4-0.5))/6.7/TD*2,abs(specs)./sum(abs(specs).^2)*6.7);
else
    plot(((0:(TD/2-1))-(TD/4-0.5))/6.7/TD*2,abs(specs)*6.7);
end
hold on
plot([0,0],ylim,'k:')
figure
hold on
plot(1:size(specd,1),specd)