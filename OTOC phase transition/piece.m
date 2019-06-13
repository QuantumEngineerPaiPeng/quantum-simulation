function y=piece(x,J,a,b,c)
y=zeros(size(x));
for p=1:length(x)
    if x(p)<J
        y(p)=c+a*(x(p)-J);
    else
        y(p)=c+b*(x(p)-J);
    end
end
end