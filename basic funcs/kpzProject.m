% 20180423
% generate the projection operator to moment k, parity P (under spatial 
% reflection), and partiy Z (under spin rotation) subspace
% k given in the unit of pi/N_atom
% P=1 or -1
% Z=1 or -1
% first colume of y gives the state, second colume gives sigma*R, third
% gives m
% see http://physics.bu.edu/~sandvik/vietri/dia.pdf pg. 22

function y=kpzProject(N_atom,k,P,Z)
y=sparse([]);

for p=1:2^N_atom
    bb=dec2bin(p-1,N_atom);
    ndown=sum(bb(:)=='1');
    if mod(ndown,2)~=(Z+1)/2
        continue
    end
    Rm=ptPeriod(p,N_atom);
    R=Rm(1);
    m=Rm(2);
    if or(R==-1,mod(k,N_atom/R)~=0) % pass if already calculated or incompatible
        continue
    end
    if or(k==0,k==N_atom/2) % sigma=+1 for k=0 or N/2
        slist=1;
    else
        slist=[1,-1]; % sigma=+-1 otherwise
    end
    for sigma=slist
        if m~=-1 % if m~=-1, sigma=+1 state and sigma=-1 state may overlap, so only take one
            if and(sigma==1,(1+sigma*P*cos(k*m*2*pi/N_atom))==0) % don't store if incompatible
                continue
            end
            if and(sigma==-1,(1-sigma*P*cos(k*m*2*pi/N_atom))~=0) % store sigma=-1 only when +1 is incompatible
                R=-1;
            end
        end
        if R>0 % store only when R>0
            col=[];
            val=[];
            for pp=0:(N_atom-1)
                coltemp=Translate(p,N_atom,pp); %T^(pp)|p>
                if sigma==1
                    valtemp=cos(k/N_atom*2*pi*pp);
                else
                    valtemp=sin(k/N_atom*2*pi*pp);
                end
                [ism,loc]=ismember(coltemp,col);
                if ~ism % create a new element if this state has not been added
                    col=[col,coltemp];
                    val=[val,valtemp];
                else % add to exisiting element
                    val(loc)=val(loc)+valtemp;
                end
                revcoltemp=bin2dec(reverse(dec2bin(coltemp-1,N_atom)))+1; %T^(pp)P|p>
                revvaltemp=valtemp*P;
                [ism,loc]=ismember(revcoltemp,col);
                if ~ism
                    col=[col,revcoltemp];
                    val=[val,revvaltemp];
                else
                    val(loc)=val(loc)+revvaltemp;
                end
            end
            if or(k==0,k==N_atom/2)
                val=val*sqrt(R)/N_atom/sqrt(2);
            else 
                val=val*sqrt(R)/N_atom;
            end
            if m~=-1
                val=val/sqrt(1+sigma*P*cos(k/N_atom*2*pi*m));
            end
            y=[y,sparse(col,ones(size(col)),val,2^N_atom,1)];
        end
    end
end
