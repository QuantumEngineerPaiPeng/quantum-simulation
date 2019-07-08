% 20190617
% Analyse 2D encoding experiment
% created by Pai Peng

n1=8; % first dimension
n2=4; % second dimension

Sall=FID(13,:);
Sall=reshape(Sall,n1*n2,[]);

Iqlist=zeros(n1,n2,size(Sall,2));

figure
ha = tight_subplot(1,size(Sall,2),[.01 .015],[.1 .01],[.01 .01]);

for p=1:size(Sall,2)
    S=Sall(:,p).';
    S=S/sqrt(S*S'); % normalize
    S=reshape(S,n1,n2);
    
    Iq=fft2(S);
    Iq=[Iq(:,(2+ceil(end/2)):end),Iq(:,1:(ceil(end/2)+1))];
    Iq=[Iq((2+ceil(end/2)):end,:);Iq(1:(ceil(end/2)+1),:)]; % shift zero freq to center
    Iqlist(:,:,p)=Iq;
    x1=(ceil(n1/2)-n1+1):(ceil(n1/2));
    x2=(ceil(n2/2)-n2+1):(ceil(n2/2));
    
    [x1,x2]=meshgrid(x1,x2);
    
    axes(ha(p))
    pcolor(x1.',x2.',abs(Iq))
    colorbar
end

set(gcf,'units','centimeters','position',[0.3,10,50,50/size(Sall,2)])