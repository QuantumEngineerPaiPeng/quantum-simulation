% 20190528
% function to change the color of a multiple lines in the same plot
% according color schemes
% Input:
% cs: color scheme. 'BW': black white, 'RB': red blue, 'rainbow' 
% Created by Pai Peng

function ppColorScheme(cs)
% cs='rainbow';

h = findobj(gca,'Type','line');
nline=length(h);

colorarray=zeros(nline,3);

if strcmp(cs,'rainbow')
    lambdalist=linspace(380,645,nline);
    for p=1:nline
        colorarray(p,:)=wavelength2RGB(lambdalist(p));
    end
else
    if strcmp(cs,'RB')
        c1=[0,0,1];
        step=[1/(nline-1),0,-1/(nline-1)];
    else
        if strcmp(cs,'BW')
            c1=[0,0,0];
            step=[0.8/(nline-1),0.8/(nline-1),0.8/(nline-1)];
        else
            error('Invalid color scheme')
        end
    end
    for p=1:nline
        colorarray(p,:)=c1+(p-1)*step;
    end
end
    
    
for p=1:nline
    h(p).Color=colorarray(p,:);
end