function outma = enlarge_pcolor(inma)
    dd1 = size(inma, 1);
    dd2 = size(inma, 2);
    outma = zeros(dd1+1, dd2+1);
    outma(1:end-1, 1:end-1) = inma;
end