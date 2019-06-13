%20180228
%simulation of the full master sequence considering the finite width of the
%pulse
%for mqc-type experiment
J=8.18e-3;
p1=1;
theta_list=[6:3:24,30:6:48];
theta_list=theta_list/360*2*pi;
q_max=3;
phi_list=(0:2*q_max-1)/(2*q_max)*2*pi;
N_atom=8;
N_rand=1;
NT=12;% number of periods
bc='p';% 'p' for periodic boundary condition; 'o' for open boundary condition
cplist=[1,1/8,1/27];
%cplist=1;
Sigma_i=IndivPauli(N_atom);
Sigma_x=sum(cat(3,Sigma_i{1}{:}),3);
Sigma_y=sum(cat(3,Sigma_i{2}{:}),3);
Sigma_z=sum(cat(3,Sigma_i{3}{:}),3);

rho_0=Sigma_y/(N_atom*2^N_atom)^0.5;
V0=Sigma_z;

%RandomM=rand(N_rand,3,N_atom);
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
        H_int=H_int+J*cplist(pp)*...
            (Sigma_i{3}{p}*Sigma_i{3}{q}-...
            (Sigma_i{1}{p}*Sigma_i{1}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
        %     H_int=H_int+(-1)*...
        %         (Sigma_i{1}{p}*Sigma_i{1}{p+1});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end

clearvars('Sigma_i')

px=expm(1i*(H_int+pi/(4*p1)*Sigma_x)*p1);% here are U^\dag
py=expm(1i*(H_int+pi/(4*p1)*Sigma_y)*p1);
pxb=expm(1i*(H_int-pi/(4*p1)*Sigma_x)*p1);
pyb=expm(1i*(H_int-pi/(4*p1)*Sigma_y)*p1);

[V,D]=eig(H_int);

%forward evolution
cycle=96;
tau=cycle/24;
u=-0.2;
v=0.2;
w=0;

t1=tau*(1-v+w)-p1/2;
t2=tau*(1-u+v)-p1;
t3=tau*(1+u-w)-p1/2;

Uf1=V*diag(exp(diag(1i*D*t1)))*V'*px*V*diag(exp(diag(1i*D*t2)))*V'*py*...
    V*diag(exp(diag(1i*D*2*t3)))*V'*py*V*diag(exp(diag(1i*D*t2)))*V'*...
    px*V*diag(exp(diag(1i*D*t1)))*V';
Uf2=Uf1;
Uf3=V*diag(exp(diag(1i*D*t1)))*V'*pxb*V*diag(exp(diag(1i*D*t2)))*V'*pyb*...
    V*diag(exp(diag(1i*D*2*t3)))*V'*pyb*V*diag(exp(diag(1i*D*t2)))*V'*...
    pxb*V*diag(exp(diag(1i*D*t1)))*V';
Uf4=Uf3;

%backward evolution
u=-0.2;
v=0;
w=0.2;

t1=tau*(1-v+w)-p1/2;
t2=tau*(1-u+v)-p1;
t3=tau*(1+u-w)-p1/2;

Ub1=V*diag(exp(diag(1i*D*t1)))*V'*px*V*diag(exp(diag(1i*D*t2)))*V'*py*...
    V*diag(exp(diag(1i*D*2*t3)))*V'*py*V*diag(exp(diag(1i*D*t2)))*V'*...
    px*V*diag(exp(diag(1i*D*t1)))*V';
Ub2=Ub1;
Ub3=V*diag(exp(diag(1i*D*t1)))*V'*pxb*V*diag(exp(diag(1i*D*t2)))*V'*pyb*...
    V*diag(exp(diag(1i*D*2*t3)))*V'*pyb*V*diag(exp(diag(1i*D*t2)))*V'*...
    pxb*V*diag(exp(diag(1i*D*t1)))*V';
Ub4=Ub3;



otoc=[];
for theta=theta_list
    
    UfT=expm(-1i*theta*Sigma_z/2)*Uf1*Uf2*Uf3*Uf4;
    Uf=UfT^(NT);
    rhotf=Uf'*rho_0*Uf;
    
    UbT=1*...
        expm(-1i*pi/4*(Sigma_z))*Ub1*Ub2*Ub3*Ub4*expm(1i*pi/4*(Sigma_z))*...
        expm(1i*theta*Sigma_z/2);
    
    Ub=UbT^(NT);
    rhotb=Ub*rho_0*Ub';
    slist=[];
    for phi=phi_list
        slist=[slist,trace(expm(-1i*phi/2*Sigma_z)*rhotf*expm(1i*phi/2*Sigma_z)*rhotb)];
    end
    
    qsquare=(0:ceil(q_max)).^2;
    qsquare=[qsquare,(q_max-1:-1:1).^2];
    Iq=fft(slist)/(2*q_max);
%     Iq=abs(Iq);
%     Iq=Iq/sum(Iq);
    otoc=[otoc,qsquare*Iq.'];
end
figure(1)
hold on
plot(theta_list/cycle/J,otoc,'LineWidth',2)