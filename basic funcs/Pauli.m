% 20180731
% Convert a single character to a 2 by 2 matrix
% Created by Pai Peng
function y=Pauli(x)
switch x
    case 'I'
        y=speye(2);
    case 'e'
        y=sparse([1,0;0,0]);
    case 'g'
        y=sparse([0,0;0,1]);
    case '+'
        y=sparse([0,1;0,0]);
    case '-'
        y=sparse([0,0;1,0]);
    otherwise
        y=sparse(eval(['sigma_',x]));
end