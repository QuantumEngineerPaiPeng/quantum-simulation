%20180301
%post processing for pulse calibration
td2=TD(2);
td1=TD(3);
Bb=abs(reshape(B,td2,td1));
bb=ones(td2,1);
for p=1:td1
    bb=bb.*Bb(:,p);
end
figure
hold on
plot(1:td2,abs(bb))

[p1,p2]=max(bb);
p2
figure
plot(1:td2,abs(Bb))