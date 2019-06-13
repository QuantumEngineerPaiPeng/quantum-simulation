%20180120
%MQC
N_atom=9;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);
bc='p';% 'p' for periodic boundary condition; 'o' for open boundary condition
cplist=[1,1/8,1/27];

J=1;
g=1.2;

T=3.77/4;
%T=0;
rho_0=Sigma_x;
q_max=3;
phi_list=(0:2*q_max-1)*2*pi/(2*q_max);

H_int=zeros(2^N_atom);
for pp=1:length(cplist)
    for p=1:N_atom-1
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else 
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
        %     H_int=H_int+(-1)*...
        %         (Sigma_i{1}{p}*Sigma_i{1}{p+1});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end

H=H_int+g*Sigma_z;

U_H=exp(1)^(-1i*H*T);
U=U_H;

rho_t=U*rho_0*U';

s_list=[];
for p=phi_list
    U_phi=exp(1)^(-1i*p/2*Sigma_z);
    rho_phi=U_phi*rho_t*U_phi';
    s_list=[s_list,trace(rho_phi*rho_t)];
end

q_list=fft(s_list)/(2*q_max);
otoc=sum(([0:q_max,q_max-1:-1:1]).^2.*q_list)/(N_atom*2^N_atom);
% figure(1)
% bar(q_list(1:q_max+1))
figure(2)
hold on
plot(g,otoc,'b+')