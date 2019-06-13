% 20180622
% Compute the overlap of 2 operator classes Tr(o1*o2)
function y=overlap(o1,o2)
if isa(o1,'double')
    y=sum(sum(o1.*(o2.')));
else
    if ~IsSameSym({o1,o2})
        error('Overlapping two operators with different length or symmetry class')
    end
    y=0;
    for p=1:length(o1.matrix)
        y=y+sum(sum(o1.matrix{p}.*(o2.matrix{p}.')));
    end
end