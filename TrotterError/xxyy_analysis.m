% 20190619
% analyse the yy and xx correlation under natural Hamiltonian trotter
% scheme
% created by Pai Peng
figure
diff=zeros(1,10);
reldiff=diff;
subplot(2,1,1)
for p=2:11
    n=round(32/p);
    plot((1:n)*p,abs(real(B(1:n,p)-B((1:n)*p,1))))
    diff(p-1)=sqrt(sum(abs(real(B(1:n,p)-B((1:n)*p,1)).^2))/n);
    reldiff(p-1)=sqrt(sum(abs(real(B(1:n,p)./B((1:n)*p,1)-1).^2))/n);
%     reldiff(p-1)=(sum(abs(real(B(1:n,p)./B((1:n)*p,1)-1))))/n;
    hold on
end
ppColorScheme('rainbow')
subplot(2,2,3)
plot(2:11,diff)
subplot(2,2,4)
plot(2:11,reldiff)