% 20180401
% find the position of an element in a 3D array
% Input: 3d logical array
% Output: 3-by-1 vector giving the index of the first appearance of True
% (first search the first dimension)
% Created by Pai Peng
function [p1,p2,p3]=find3d(A)
s=size(A);
for p1=1:s(1)
    for p2=1:s(2)
        for p3=1:s(3)
            if A(p1,p2,p3)
                return
            end
        end
    end
end
