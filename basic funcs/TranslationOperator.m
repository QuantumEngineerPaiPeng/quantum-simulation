% 20180423
% generate translation operator for given number of atoms and translate
% distance
function T=TranslationOperator(N_atom,distance)
row=[];
for p=0:2^N_atom-1
    pp=dec2bin(p,N_atom);
    ppp=[pp(distance+1:end),pp(1:distance)];
    row=[row,bin2dec(ppp)+1];
end
T=sparse(row,1:2^N_atom,1);
end