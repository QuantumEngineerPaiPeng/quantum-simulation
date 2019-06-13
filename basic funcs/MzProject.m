% 20180913
% For a given L, classify all the 2^L basis according to Mz
% Input: L
% Output: L+1 cell, nth element contains the label of state with Mz=-L-2+2n
% Created by Pai Peng
function y=MzProject(L)
if isunix
    filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/MzSym_%d.mat',L);
else
    filename=['C:\Users\Pai\Dropbox (MIT)\grad\research\codes\basic funcs\data\',sprintf('MzSym_%d.mat',L)];
end
if exist(filename,'file')
    load(filename,'y')
else
    y=cell(L+1,1);
    for p=1:2^L
        ps=dec2base(p-1,2,L);
        ps(ps=='1')=[];
        d=length(ps);
        y{d+1}=[y{d+1},p];
    end
    save(filename,'y')
end