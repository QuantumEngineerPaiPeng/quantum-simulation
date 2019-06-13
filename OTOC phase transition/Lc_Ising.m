%20171016
%Compute correlation length of transverse Ising model, with density matrix
%analytically obtained
N_spin=50;
U=zeros(N_spin,N_spin);
V=U;
F=U;
J=1;
g=0.5;
%t=0.16*12;

min_t=0;
max_t=10;
N_t=100;
step_t=(max_t-min_t)/N_t;
x=min_t:step_t:max_t;
L_c=zeros(1,N_t+1);
L_c1=L_c;%L_c1 is when (r+1)*(N-r) is replaced by (r+1)*N
plist=-(N_spin-1)/2:1:(N_spin-1)/2;
plist=plist*2*pi/N_spin;

count=1;
for t=min_t:step_t:max_t
    p=1;
    
    for p=1:N_spin
        w_p=(J^2+g^2+2*J*g*cos(plist(p)))^0.5;
        sin_theta_p=J*sin(plist(p))/w_p;
        cos_theta_p=(J*cos(plist(p))+g)/w_p;
        U(p,p)=cos(w_p*t)-1i*cos_theta_p*sin(w_p*t);
        V(p,N_spin+1-p)=-sin_theta_p*sin(w_p*t);
    end
    
    for row=1:N_spin
        for column=1:N_spin
            F(row,column)=exp(1i*column*plist(row))/N_spin^0.5;
        end
    end
    
    u=F'*U*F;
    v=transpose(F)*V*F;
    
    rho_0=u'*u-v'*v;
    rho_2=-u'*v-v*u';
    rho_0=rho_0/(2*N_spin)^0.5;
    rho_2=rho_2/(2*N_spin)^0.5;
    
    for p=-N_spin:N_spin
        L_c(count)=L_c(count)+(abs(p)+1)*sum(abs(diag(rho_0,p)).^2);
        L_c(count)=L_c(count)+(abs(p)+1)*sum(abs(diag(rho_2,p)).^2);
    end
    L_c1(count)=N_spin*abs(rho_0(1,1))^2;
    for p=2:N_spin
        L_c1(count)=L_c1(count)+(N_spin-p+1)*(2*p*abs(rho_0(1,p))^2+2*p*abs(rho_2(1,p))^2);
    end
    count=count+1;
end
figure(1)
hold on
L_c=L_c*2;
L_c1=L_c1*2;
plot(x,L_c)%,x,left_limit(J,g,x))