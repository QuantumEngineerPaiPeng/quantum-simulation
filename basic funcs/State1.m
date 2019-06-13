% 20180820
% function to generate single spin-1/2 pure state
% Created by Pai Peng
function y=State1(dir)
if isa(dir,'char')
    switch dir
        case 'z'
            y=[1;0];
        case '-z'
            y=[0;1];
        case 'x'
            y=[1;1]/sqrt(2);
        case '-x'
            y=[1;-1]/sqrt(2);
        case 'y'
            y=[1;1i]/sqrt(2);
        case '-y'
            y=[1;-1i]/sqrt(2);
        otherwise
            error('argument invalid')
    end
else
    if isa(dir,'double')
        y=[cos(dir(1)/2/180*pi);sin(dir(1)/2/180*pi)*exp(1i*(dir(2)/180*pi))];
    else
        error('argument invalid')
    end
end