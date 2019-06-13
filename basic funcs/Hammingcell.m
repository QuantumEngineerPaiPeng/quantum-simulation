% 20180830
% For a given L, classify all d^L prodect operators according to their
% Hamming weight, with d being the number of basis. d=4 when all Pauli
% operators are involved, d=2 when only z is involved.
% Output: cell. nth element is an array that contains all the Hamming
% weight=n-1. In the array, each number n represents a product operator that
% is given by n-1 in base d, e.g 1 stands for 0000, 6 stands for 0011 in
% d=4, L=4.
% The cell is saved in basic funcs -> data.
% Created by Pai Peng
function y1=Hammingcell(L,d)
if isunix
    filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/Hammingcell%d_%d.mat',L,d);
    if ~exist(filename,'file')
        y=cell(1,L+1);
        for p=1:(d^L)
            ps=dec2base(p-1,d);
            ps(ps=='0')=[]; % delete all '0' terms
            y{length(ps)+1}=[y{length(ps)+1},p];
        end
        save(filename,'y')
    else
        load(filename)
    end
else
    filename=['C:\Users\Pai\Dropbox (MIT)\grad\research\codes\basic funcs\data\',sprintf('Hammingcell%d_%d.mat',L,d)];
    if ~exist(filename,'file')
        y=cell(1,L+1);
        for p=1:(d^L)
            ps=dec2base(p-1,d);
            ps(ps=='0')=[]; % delete all '0' terms
            y{length(ps)+1}=[y{length(ps)+1},p];
        end
        save(filename,'y')
    else
        load(filename)
    end
end
y1=y;
end