% 20180905
% For a given 2^L-by-2^L matirx, each matrix element correspond to a
% product of |e><e|, |g><g|, |e><g|, |g><e|. Classify all matrix elements
% according to how many |e><g| and |g><e| are involved,
% Output: cell. nth element is an array that contains all the Hamming
% weight=n-1. In the array, each number n represents a product operator that
% is given by n-1 in base d, e.g 1 stands for 0000, 6 stands for 0011 in
% L=4, where 0 represents |e><e|, 1: |g><g|, 2: |e><g|, 3: |g><e|. 2 and 3
% are off diagonal.
% The cell is saved in basic funcs -> data.
% Created by Pai Peng
function y1=OffDiagcell(L)
filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/OffDiagcell%d_%d.mat',L,d);
if ~exist(filename,'file')
    y=cell(1,L+1);
    for p=1:(d^L)
        ps=dec2base(p-1,d);
        ps(ps=='0')=[]; % delete all '0' terms
        ps(ps=='1')=[]; % delete all '1' terms
        y{length(ps)+1}=[y{length(ps)+1},p];
    end
    save(filename,'y')
else
    load(filename)
end
y1=y;
end