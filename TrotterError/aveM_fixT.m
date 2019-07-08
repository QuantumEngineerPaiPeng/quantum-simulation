% 20190426
% Calculate averaged M as a function of rotation angle from experimental
% data. Averaged upto the total experiment time of the experiment with
% smallest angle.
% The experimental data in 'B' is taken for the same number of periods for
% different rotation angle, so we need to truncate the angle to get the
% average.
% Created by Pai Peng
philist=15:15:150; % rotation angle of the entire data set store in 'B'
% philist=([2,4,8,16,24,32,40,48,56,64,72,80]+2.95)*4*8.18e-3/pi*180;
% philist=([2,4,8,16,32,64]+2.95)*4*8.18e-3/pi*180;
explist=1:10; % subset of the data to analysis.  
% Use the entire set may not be good if the max(phi)/min(phi) is too large,
% bcs the largest phi is only averaged over a small number of times.
philist1=philist(explist);
B1=abs(B(:,explist));
Nlist=round(size(B,1)./(philist/min(philist))); % # of periods to be averaged

Bmean1=zeros(length(explist),1);
Bmean2=zeros(length(explist),1);
Bt=zeros(length(explist),1);
for p=1:length(explist)
    Bmean1(p)=mean(B1(1:Nlist(p),p));
    Bmean2(p)=mean(B1(round(Nlist(p)/2):Nlist(p),p));
    Bt(p)=B1(Nlist(p),p);
end

figure
plot(philist1,[Bmean1,Bmean2,Bt])
xlabel('rotation angle $h\tau$ (degree)','Interpreter','latex')
legend({'ave. 0 to T','ave. T/2 to T','at T'})
hold on
ppStyle(30,2,10)