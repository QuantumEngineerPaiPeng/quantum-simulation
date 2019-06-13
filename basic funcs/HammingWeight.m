% 20180817
% Recursion function to calculate the Hamming weight of a 2^L by 2^L matrix
% If the input matrix is 2 by 2, then returns 2 numbers giving the
% 0-spin weight and 1-spin weight.
% If the input matrix is in larger dimension, say 2^L by 2^L, decompose M
% to kron(s0,M0), kron(s1,M1), kron(s2,M2), kron(s3,M3) and return 
% [HammingWeight(M0),0]+[0,HammingWeight(M1)+HammingWeight(M2)+HammingWeight(M3)]
% Created by Pai Peng
function y=HammingWeight(M)
if isscalar(M)
    y=M*M';
    return
end
L=size(M,1)/2;
if L~=size(M,2)/2
    error('Input must be a square matrix')
end
M11=M(1:L,1:L);
M12=M(1:L,(L+1):end);
M21=M((L+1):end,1:L);
M22=M((L+1):end,(L+1):end);

y=[HammingWeight((M11+M22)/2),0]+...
    [0,HammingWeight((M12+M21)/2)+HammingWeight((M12-M21)/(-2*1i))+HammingWeight(((M11-M22)/2))];