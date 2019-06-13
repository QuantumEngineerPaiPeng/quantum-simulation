% 20180423
% translate a produce space basis specified by an integer in decimal by
% distance
function y=Translate(p,N_atom,distance)
    pp=dec2bin(p-1,N_atom);
    ppp=[pp(distance+1:end),pp(1:distance)];
    y=bin2dec(ppp)+1;
end