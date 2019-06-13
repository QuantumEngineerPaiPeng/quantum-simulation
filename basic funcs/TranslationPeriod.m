% 20180423
% generate the tranlation period of base vectors in product space
function y=TranslationPeriod(N_atom)
y=zeros(2^N_atom,1);
for p=1:2^N_atom
    for distance=1:(N_atom)
        pp=Translate(p,N_atom,distance);
        if pp<p
            y(p)=-1;
            break
        end
        if pp==p
            y(p)=distance;
            break
        end
    end
end
end