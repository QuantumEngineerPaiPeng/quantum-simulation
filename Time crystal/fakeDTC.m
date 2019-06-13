% 20180417
% fake DTC data to understand CF responce to decay rate
Nlist=8:8:200;
cf=zeros(size(Nlist));
count=1;
for N=Nlist
gamma=0.02;
B1=mod(0:(N-1),2)-0.5;
B2=exp(gamma*(1:N));
B=B2.*B1;
F=fft(B);
cf(count)=abs(F(length(B)/2+1))^2/sum(abs(F).^2);
count=count+1;
end
figure(1)
plot(Nlist,cf)