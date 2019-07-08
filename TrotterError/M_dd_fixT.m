% 20190429
% Get the magnetization, as a function of rotation angle from experimental
% data. 'expt3' 41-49
% Created by Pai Peng
philist=15:15:165; % rotation angle 
explist=[381:391]; % experiment number
Mlist1=zeros(size(explist)); % at a given t
Mlist2=zeros(size(explist)); % averaged over the later half of the experiments
figure
for p=1:length(explist)
    data_analysis(explist(p),'expt3')
    plot(1:length(B),real(B))
    hold on
    Mlist1(p)=real(B(end));
    Mlist2(p)=real(B(round(end/4)));
%     Mlist2(p)=mean(real(B(round(end/2):end)));
end
xlabel('# of periods')

figure
hold on
plot(philist,Mlist1)
xlabel('h\tau')
ppStyle(30,2,10)
figure
hold on
plot(philist,Mlist2)
xlabel('h\tau')
ppStyle(30,2,10)