% 20180316
% Data analysis for 2D DTC experiment
bb=permute(B,[1,3,2]);

sh=size(bb);
bb=reshape(bb,sh(1)*sh(2),sh(3));

row=5;
col=4;

figure(1)
for p=1:row*col
    subplot(row,col,p)
    plot(real(bb(:,p)),'LineWidth',1)
end