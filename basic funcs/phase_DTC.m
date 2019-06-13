%20180307
%generate pulse phase list for DTC with pi pulse implemented as phase shift
L=30;
epsilon=42;
sign='+';% '-' for pi-epsilon '+' for pi+epsilon
phlist=[];

if sign=='-'
    name=['ph0';'ph2';'ph1';'ph3'];
    iph=[0,180,90,270];
    for p=1:4
        line=[name(p,:),' = (360)'];
        for pp=0:L-1
            line=[line,' ',num2str(mod(iph(p)+pp*(180-epsilon),360))];
        end
        phlist=[phlist,line,newline];
    end
else
    name=['ph0';'ph2';'ph1';'ph3'];
    iph=[0,180,90,270];
    for p=1:4
        line=[name(p,:),' = (360)'];
        for pp=0:L-1
            line=[line,' ',num2str(mod(iph(p)-pp*(180-epsilon),360))];
        end
        phlist=[phlist,line,newline];
    end
end

disp(phlist)