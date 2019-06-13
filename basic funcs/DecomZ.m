% 20180829
% Recursion function to decompose a 2^L-by-2^L diagonal matrix to the
% products of sz
% If the input matrix is 1-by-1, then returns the scalar.
% If the input matrix is in larger dimension, say 2^L by 2^L, decompose M
% to kron(s0,M0), kron(sz,Mz) then compute
% DecomZ(M0) and DecomZ(Mz), and append the resulting array
% together.
% Created by Pai Peng
function y=DecomZ(M)
if isscalar(M)
    y=M;
    return
end
if size(M,1)==1
    L=size(M,2)/2;
else
    if size(M,2)==1
        L=size(M,1)/2;
    else
        error('Input must be a vector')
    end
end

M1=M(1:L);
M2=M((L+1):end);
y=[DecomZ((M1+M2)/2),DecomZ(((M1-M2)/2))];