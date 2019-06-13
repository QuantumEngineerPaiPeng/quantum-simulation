% 20180906
% Recursion function to decompose a 2^L by 2^L matrix to the products of
% {|e><e|,|g><g|,|e><g|,|g><e|}
% If the input matrix is 1 by 1, then returns the scalar.
% If the input matrix is in larger dimension, say 2^L by 2^L, decompose M
% to kron(|e><e|,M0), kron(|g><g|,M1), kron(|e><g|,M2), kron(|g><e|,M3) and compute
% DecomEG(M0) to EG(M3), and append the resulting array
% together.
% Created by Pai Peng
function y=DecomEG(M)
if isscalar(M)
    y=M;
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
y=[DecomEG(M11),DecomEG(M22),DecomEG(M12),DecomEG(M21)];