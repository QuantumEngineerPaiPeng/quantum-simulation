% 20180423
% generate the projection operator to moment k subspace
function y=kProject(N_atom,k)
R=TranslationPeriod(N_atom);
y=sparse([]);
for p=1:2^N_atom
    if or(R(p)==-1,mod(k,N_atom/R(p))~=0)
        continue
    end
    col=p;
    val=1;
    for pp=1:(N_atom-1)
        coltemp=Translate(p,N_atom,pp);
        [ism,loc]=ismember(coltemp,col);
        if ~ism
            col=[col,coltemp];
            val=[val,exp(1i*k/N_atom*2*pi*pp)];
        else
            val(loc)=val(loc)+exp(1i*k/N_atom*2*pi*pp);
        end
    end
    val=val*sqrt(R(p))/N_atom;
    y=[y,sparse(col,ones(size(col)),val,2^N_atom,1)];
end