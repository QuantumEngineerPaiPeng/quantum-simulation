% 20180417
% extract crystalline fraction from experimental data
% data stored in variable B, first dim is time, second is perturbation or
% period, specifie by varname
Bb=real(squeeze(B));
varname='perturbation';
perf_pi=1.91;
Xlist=(1.71:0.05:2.11)-perf_pi;
%seg={1:32,33:64,65:96,97:128};
% seg={1:64,65:128};
% seg={1:32,33:64};
seg={1:64};
cf=zeros(size(B,2),length(seg));
for p=1:size(B,2)
    for pp=1:length(seg)
        M=B(seg{pp},p);
        F=fft(M);
        cf(p,pp)=abs(F(length(M)/2+1))^2/sum(abs(F).^2);
    end
end

figure(1)
hold on
plot(Xlist,cf)
xlabel(varname)
ylabel('crystalline fraction')
ppStyle(30,2)