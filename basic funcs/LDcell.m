% 20180910
% For a given L, classify all d^L prodect operators according to their
% correlation distance to certain spin s, with d being the number of basis.
% d=4 when all Pauli operators are involved, d=2 when only z is involved.
% Output: cell. nth element is an array that contains all the correlation
% distance=n-1 terms. In the array, each number n represents a product operator that
% is given by n-1 in base d, e.g 1 stands for 0000, 6 stands for 0011 in
% d=4, L=4.
% The cell is saved in basic funcs -> data.
% Created by Pai Peng
function y=LDcell(L,d,s)
if isunix
    filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/data/LDcell%d_%d_%d.mat',L,d,s);
    if ~exist(filename,'file')
        y=cell(1,L);
        for p=2:d^L
            ps=dec2base(p-1,d,L);
            for q=1:L
                if ps(q)~='0'
                    s1=q;
                    break
                end
            end
            
            for q=L:-1:1
                if ps(q)~='0'
                    s2=q;
                    break
                end
            end
            
            D=max(abs([s1-s,s2-s]));
            y{D+1}=[y{D+1},p];
        end
        save(filename,'y')
    else
        load(filename)
    end
else
    filename=['C:\Users\Pai\Dropbox (MIT)\grad\research\codes\basic funcs\data\',sprintf('LDcell%d_%d_%d.mat',L,d,s)];
    if ~exist(filename,'file')
        y=cell(1,L);
        for p=2:d^L
            ps=dec2base(p-1,d,L);
            for q=1:L
                if ps(q)~='0'
                    s1=q;
                    break
                end
            end
            
            for q=L:-1:1
                if ps(q)~='0'
                    s2=q;
                    break
                end
            end
            
            D=max(abs([s1-s,s2-s]));
            y{D+1}=[y{D+1},p];
        end
        save(filename,'y')
    else
        load(filename)
    end
end

