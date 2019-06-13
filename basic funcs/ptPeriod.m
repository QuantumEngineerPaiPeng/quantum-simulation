% 20180423
% generate the tranlation period R of base vectors in product space
% and the reflection-translation (T^m*P) period m
function y=ptPeriod(p,N_atom)
    m=-1;
    for distance=1:(N_atom)
        pp=Translate(p,N_atom,distance);
        if pp<p
            R=-1;
            m=-1;
            break
        end
        if pp==p
            R=distance;
            m=-1;
            break
        end
    end
    revp=bin2dec(reverse(dec2bin(p-1,N_atom)))+1;
    pp=revp;
    for distance=0:R-1
        if pp<p
            R=-1;
            m=-1;
            break
        else
            if pp==p
                m=distance;
                break;
            end
        end
        pp=Translate(pp,N_atom,1);
    end
    y=[R,m];
end