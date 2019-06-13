% 20180830
% For a given cell tau, get the overlap between every two elements
% Created by Pai Peng
function y=GetOverlapMat(tau)
y=zeros(length(tau),length(tau));
for p=1:length(tau)
    for q=p:length(tau)
        y(p,q)=overlap(tau{p},tau{q});
        y(q,p)=y(p,q);
    end
end