% 20180321
% Data analysis for window DTC experiment
Bb=abs(B)./abs(B(1,:));
Bb=Bb';
sh=size(Bb);
figure(1)
hold on
plot(1:sh(2),Bb(:,:),'LineWidth',2)