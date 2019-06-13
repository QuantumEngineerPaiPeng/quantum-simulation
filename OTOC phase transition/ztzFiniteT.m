%20180213
%compute the ztz averaged OTOC at any temperature
ome=@(p,g) sqrt(1-2*g*cos(p)+g^2);
fun=@(p,g,beta) (16/pi*sin(p).^2./(ome(p,g).^2)-12/pi*sin(p).^4./(ome(p,g).^4))./(1+(cosh(beta*ome(p,g))).^(-1));
g_list=[0.1:0.05:3];
%beta_list=logspace(-2,2,10);
beta_list=0;
C=zeros(length(g_list),length(beta_list));

row=1;
for g=g_list
    col=1;
    for beta=beta_list
        C(row,col)=integral(@(x)fun(x,g,beta),0,2*pi);
        col=col+1;
    end
    row=row+1;
end

figure(1)
plot(g_list,C)