% 20180829
% function to reshape a 1-by-2^L array to 2-by-2-by-2...by-2 tensor
function y=reshape2L(x)
L=round(log2(length(x)));
y=eval(['reshape(x,2',repmat(',2',1,L-1),');']); 
end