% 20180711
% Compute Tr|M|^2, with M given by an OperatorClass
% Created by Pai Peng
function y=TwoNorm(rho)
y=0;
for p=1:length(rho.matrix)
    y=y+full(sum(sum(abs(rho.matrix{p}).^2)));
end