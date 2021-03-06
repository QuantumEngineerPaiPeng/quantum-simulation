%20170917
%compute c(t)=<[Y(t),X]> with different epsilon with fixed interaction J
%find the peak eps_M(J). Then vary J find the explicit function eps_M(J).
%See how the disorder affect this function.
N_atom=8;
N_rand=10;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

Disorder=[0,25,0]; %Random Lamor freq. is uniformly discributed in [-0.5*Disorder, 0.5*Disorder]
tau=1;

t=12;

rho_0=Sigma_y;

min_J=0.02;
N_J=10;
max_J=0.152;
step_J=(max_J-min_J)/(N_J-1);

min_eps=0.03;
N_eps=50;
max_eps=0.1;
step_eps=(max_eps-min_eps)/(N_eps-1);

%RandomM=rand(N_rand,3,N_atom);
H_int=zeros(2^N_atom);
C=zeros(N_J,N_eps,N_rand);

for p=1:N_atom-1
    H_int=H_int+(-1)*...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
%         H_int=H_int+(-2)*...
%             (Sigma_i{2}{p}*Sigma_i{2}{p+1});
end

for pp=1:N_rand
    count_J=1;
    H_rand=zeros(2^N_atom);
    for p=1:N_atom
        H_rand=H_rand+Disorder(1)*(rand-0.5)*Sigma_i{1}{p}+...
            Disorder(2)*(rand-0.5)*Sigma_i{2}{p}+...
            Disorder(3)*(rand-0.5)*Sigma_i{3}{p};
    end
    for J=min_J:step_J:max_J
        g=0.152-J;%strength of the disorder
%         g=0;
        count_eps=1;
        for epsilon=min_eps:step_eps:max_eps;
            
            U_p=exp(1)^(-1i*Sigma_x*pi/2*(1-epsilon));
%             U_p_inv=exp(1)^(1i*Sigma_x*pi/2*(1-epsilon));
            H=J*H_int+g*H_rand;
            U_H=exp(1)^(-1i*H*tau);
%             U_H_inv=exp(1)^(1i*H*tau);
            U=U_p*U_H;
%             U_inv=U_H_inv*U_p_inv;
            rhot=(U)^t*rho_0*(U)^-t;
            C(count_J,count_eps,pp)=-trace((rhot*Sigma_x-Sigma_x*rhot)^2);
            count_eps=count_eps+1;
        end
        count_J=count_J+1;
    end
    pp
end

C_eps_J=real(squeeze(sum(C,3)))'/N_rand;
X=min_eps:step_eps:max_eps;

if N_rand==1
    figure(3)
    plot(X,C_eps_J)
else
    stdev=std(real(C),1,3)'/N_rand^0.5;
    figure(3)
    hold on
    for p=1:N_J
        errorbar(X,C_eps_J(:,p),stdev(:,p))
    end
end
[C_M,Ind]=max(C_eps_J);
eps_M=X(Ind);
Y=min_J:step_J:max_J;
figure(4)
plot(eps_M,Y)