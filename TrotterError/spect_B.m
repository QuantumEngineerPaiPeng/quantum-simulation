% 20190425
% get spectrum from variable 'B'
normalize=false;
dtlist=[([2,4,8,16,32]+2.95)*2,64+2.95];
% dtlist=(32+2.95)*2;
B1=imag(B); % construct complex data
spec=fft(B1,size(B1,1),1); 
specs=[spec((end/2+1):end,:);spec(1:end/2,:)];
specd=abs(spec((end/2+1):end,:))-abs(spec(end/2:-1:1,:));
% specd=
figure(5)
for p=1:size(B1,2)
    if normalize
        plot(((0:(size(B1,1)-1))-(size(B1,1)-1)/2)/dtlist(p)/size(B1,1),abs(specs(:,p))./sum(abs(specs(:,p)).^2)*dtlist(p));
    else
        plot(((0:(size(B1,1)-1))-(size(B1,1)-1)/2)/dtlist(p)/size(B1,1),abs(specs(:,p))*dtlist(p));
    end
    hold on
end
