%20190515
%generate pulse phase list for master sequence where the phases are
%inserted every 4 pulses
N=24; % unit=360/N
N_period=0; % starting number of periods
L=24; % length of the phase list
step=1; % phase shift, in unit of 360/N degree
phlist=[];
for seq=['i','d'] %'d' for decreasing order, 'i' for increasing
    
    if seq=='d'
        iph=[N/4,3*N/4,N/2,0]+(N_period)*step;
        name=['ph0';'ph2';'ph3';'ph1'];
    else
        iph=[0,N/2,N/4,3*N/4];
        name=['ph10';'ph12';'ph11';'ph13'];
    end
    
    
    for p=1:4
        if seq=='d'
            A=0:-step:0-(step*(L-1));
        else
            A=0:step:0+(step*(L-1));
        end
        
%         A=kron(A,[1,1]);
%         A=A(2:(end-1));
        A=A.*cos(pi*(1:L));
        
        if mod(p,2)==1
            A=A+iph(p);
        else 
            A=-[A(2:end),A(1)]+iph(p);
        end
        
        line=[name(p,:),' = ','(',num2str(N),')'];
        for a = A
            line=[line,' ',num2str(mod(a,N))];
        end
        phlist=[phlist,newline,line];
    end
    if seq=='d'
        phlist=[phlist,newline];
        xphase=-step; % xphase=-step/2 for 2nd Trotter and -step for Trotter
        phlist=[phlist,newline,'ph6 = (120) ',num2str(mod(xphase,120)),' '...
            ,num2str(mod(xphase,120)),' ',num2str(mod(xphase,120)),' ',num2str(mod(xphase,120)),' '...
            ,num2str(mod(xphase+60,120)),' ',num2str(mod(xphase+60,120)),' '...
            ,num2str(mod(xphase+60,120)),' ',num2str(mod(xphase+60,120)),' '];
        
        phlist=[phlist,newline,'ph7 = (120) ',num2str(mod(xphase+60,120)),' ',num2str(mod(xphase+60,120)),' '...
            ,num2str(mod(xphase+60,120)),' ',num2str(mod(xphase+60,120)),' ',num2str(mod(xphase,120)),' '...
            ,num2str(mod(xphase,120)),' ',num2str(mod(xphase,120)),' ',num2str(mod(xphase,120)),' '...
            ];
    end
    phlist=[phlist,newline];
end
disp(phlist)