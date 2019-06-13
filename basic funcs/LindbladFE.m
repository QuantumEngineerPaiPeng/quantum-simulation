function rho=LindbladFE(rho0,H,gamma,L,Dt,N)
rho=rho0;
for p=1:N
    rho=rho+Dt*(-1i*(H*rho-rho*H)+gamma*L*rho*L'-0.5*gamma*(L'*L)*rho-0.5*gamma*rho*(L'*L));
end