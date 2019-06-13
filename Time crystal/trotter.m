% test Trotter expansion
N_atom=6;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

H1=zeros(2^N_atom);
for p=1:2:N_atom-1
    H1=H1+(-2)*...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
end

H2=zeros(2^N_atom);
for p=2:2:N_atom-1
    H2=H2+(-2)*...
        (Sigma_i{2}{p}*Sigma_i{2}{p+1}-...
        (Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1})/2);
end

H3=0.5*Sigma_z;
T=40;
N_list=logspace(0,10,50);
err1=[];
M_err1=[];
M_err2=[];
err2=[];
C_err1=[];
C_err2=[];
rho_err1=[];
rho_err2=[];
exact=exp(1)^(1i*(H1+H2+H3)*T);
rho=Sigma_z;
rhot_exact=exact*rho*exact';
M_exact=trace(rho*rhot_exact);
C_exact=-trace((rhot_exact*rho-rho*rhot_exact)^2);
for N=N_list
    t=T/N;
    first=(exp(1)^(1i*H1*t)*exp(1)^(1i*H2*t)*exp(1)^(1i*H3*t))^(N);
    second=(exp(1)^(1i*H1*t/2)*exp(1)^(1i*H2*t/2)*exp(1)^(1i*H3*t)*exp(1)^(1i*H2*t/2)*exp(1)^(1i*H1*t/2))^N;
    rhot_first=first*rho*first';
    rhot_second=second*rho*second';
    M_first=trace(rho*rhot_first);
    M_second=trace(rho*rhot_second);
    C_first=-trace((rhot_first*rho-rho*rhot_first)^2);
    C_second=-trace((rhot_second*rho-rho*rhot_second)^2);
    M_err1=[M_err1,abs(M_first-M_exact)];
    M_err2=[M_err2,abs(M_second-M_exact)];
    C_err1=[C_err1,abs(C_first-C_exact)];
    C_err2=[C_err2,abs(C_second-C_exact)];
    err1=[err1,sqrt(trace((exact-first)*(exact-first)'))];
    err2=[err2,sqrt(trace((exact-second)*(exact-second)'))];
    rho_err1=[rho_err1,sqrt(trace((rhot_exact-rhot_first)*(rhot_exact-rhot_first)'))];
    rho_err2=[rho_err2,sqrt(trace((rhot_exact-rhot_second)*(rhot_exact-rhot_second)'))];
end
figure(1)
loglog(N_list,C_err1,N_list,C_err2)
figure(2)
loglog(N_list,rho_err1,N_list,rho_err2)