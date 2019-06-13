% 20180820
% Generate the state or density matrix of a product state
% Input: dir, cell, containing individual spin state
% opt: 'state' or 'ope', optional, default 'ope'
% Created by Pai Peng
function y=StateN(dir,opt)
v=1;
for p=1:length(dir)
    v=kron(v,State1(dir{p}));
end
if nargin==2 && strcmp(opt,'state')
    y=v;
else
    y=OperatorClass(length(dir));
    y.matrix={v*v'};
end