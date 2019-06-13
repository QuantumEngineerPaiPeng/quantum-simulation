% 20180620
% Judge whether or not the symmetry class of the operators are the same
% Argument is given by a cell containing the operators.
% created by Pai Peng
function y=IsSameSym(opecell)
if length(opecell)<2
    error('Less than 2 operators are given for comparison')
end
L=opecell{1}.L;
sym=opecell{1}.sym;
for p=2:length(opecell)
    if opecell{p}.L~=L|| (~strcmp(opecell{p}.sym,sym)&& ~(isempty(opecell{p}.sym) && isempty(sym)))
        y=false;
        return
    end
end
y=true;
end