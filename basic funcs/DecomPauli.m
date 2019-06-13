% 20180730
% Recursion function to decompose a 2^L by 2^L matrix to the products of Pauli operators
% If the input matrix is 1 by 1, then returns the scalar.
% If the input matrix is in larger dimension, say 2^L by 2^L, decompose M
% to kron(s0,M0), kron(s1,M1), kron(s2,M2), kron(s3,M3) and compute
% DecomPauli(M0) to DecomPauli(M3), and append the resulting array
% together.
% Created by Pai Peng
function y=DecomPauli(M)
if isscalar(M)
    y=M;
    %y=[(M(1,1)+M(2,2))/2,(M(1,2)+M(2,1))/2,(M(1,2)-M(2,1))/(-2*1i),(M(1,1)-M(2,2))/2];
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
y=[DecomPauli((M11+M22)/2),DecomPauli((M12+M21)/2),...
    DecomPauli((M12-M21)/(-2*1i)),DecomPauli(((M11-M22)/2))];
% y=[DecomPauli((M11+M22)/2),DecomPauli(M12),...
%     DecomPauli(M21),DecomPauli(((M11-M22)/2))];