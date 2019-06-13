% 20181005
% For a given L, classify all d^L prodect operators according to the number
% of (not 0 nor m) operators, with d being the number of basis. 
% Output: cell. nth element is an array that contains all the number=n-1
% terms. In the array, each number n represents a product operator that 
% is given by n-1 in base d, e.g 1 stands for 0000, 6 stands for 0011 in
% d=4, L=4.
% The cell is saved in basic funcs -> data.
% Created by Pai Peng
function y1=notmcell(L,d,m)
if isunix
    filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/notmcell%d_%d_%d.mat',L,d,m);
else
    filename=['C:\Users\Pai\Dropbox (MIT)\grad\research\codes\basic funcs\data\',sprintf('notmcell%d_%d_%d.mat',L,d,m)];
end
if ~exist(filename,'file')
    y=cell(1,L+1);
    for p=1:(d^L)
        ps=dec2base(p-1,d);
        ps(ps=='0')=[];
        ps(ps==num2str(m))=[];
        y{length(ps)+1}=[y{length(ps)+1},p];
    end
    save(filename,'y')
else
    load(filename,'y')
end

y1=y;