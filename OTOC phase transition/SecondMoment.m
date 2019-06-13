% 20190322
% Code to calculate the second momentum of Iq by assuming equal probability
% of all operators
% Created by Pai Peng

Llist=5:50; % number of non-trivial Puali operators
fm=zeros(length(Llist),1);
sm=fm;
for pL=1:length(Llist)
    L=Llist(pL);
    
    % a1 is the number of sigma_+, a1-q is the number of sigma_-, L-2*a1+q is
    % the number of sigma_z.
    Pqa1=@(q,a1) nchoosek(L,a1)*(1/3)^a1*nchoosek(L-a1,a1-q)*(1/3)^(a1-q)*(1/3)^(L-2*a1+q); %P(Q=q and A1=a1)
    
    Pq=zeros(L+1,1); % vector that store P(Q=q)=\sum_{a1=q}^floor((L+q)/2) P(Q=q and A1=a1), with 0<=q<=L
    
    for q=0:L
        for a1=q:floor((L+q)/2)
            Pq(q+1)=Pq(q+1)+Pqa1(q,a1);
        end
    end
    
    fm(pL)=sum(Pq.*(0:L).'); % first moment of Pq (only the q>0 part)
    sm(pL)=2*sum(Pq.*(0:L).'.^2); % second mement of Pq
end

figure
plot(Llist,[fm,sm])
ppStyle(20,2,10)