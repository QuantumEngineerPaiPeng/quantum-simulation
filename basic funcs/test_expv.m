% 20190706
% test the perfomance of expv
% created by Pai Peng
tic
N_atom=13;
H_int=OperatorClass(N_atom,[0.5,0.5,-1],1,'p',1);
v=ones(2^N_atom,1);
v=v/sqrt(v'*v);
rho_0=OperatorClass(N_atom);
rho_0.matrix={v*v'};
tlist=logspace(-1,2,20);
% tollist=logspace(-10,-5,20);
mlist=5:5:100;
timelist=zeros(1,length(tlist));
overlap1=timelist;

w=v;
% for p=1:length(tlist)
%     t1=cputime;
%     w = expv( tlist(p), -1i*H_int.matrix{1}, v, 1e-5, 30 );
%     t2=cputime-t1;
%     timelist(p)=t2;
%     overlap1(p)=abs(v'*w)^2;
% end

% for p=1:length(tollist)
%     t1=cputime;
%     w = expv( 20, -1i*H_int.matrix{1}, v, tollist(p), 30 );
%     t2=cputime-t1;
%     timelist(p)=t2;
%     overlap1(p)=abs(v'*w)^2;
% end

for p=1:length(mlist)
    t1=cputime;
    w = expv( 100, -1i*H_int.matrix{1}, v, 1e-7, mlist(p) );
    t2=cputime-t1;
    timelist(p)=t2;
    overlap1(p)=abs(v'*w)^2;
end

% t1=cputime;
overlap2=TPC(rho_0,H_int,rho_0,100);
% t2=cputime-t1;

figure(2)
hold on
% loglog(tlist,timelist)
% figure
plot(mlist,abs(overlap1-overlap2))
toc