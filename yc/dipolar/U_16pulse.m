function U_real = U_16pulse( ...,
    N_atom, bc, sym, Jz, u_list, tau, hz ...,
)
% 2019 pi
% return Floquet U of 16 pulse sequence (dip_y)
% [t1 x t2 y 2t1 y t2 x, 2t1 x t2 y 2t1 y t2 x, 2t1 -x t2
% -y 2t1 -y t2 -x, 2t1 -x t2 -y 2t1 -y t2 -x t1],
% x:=e^{-ipi X}, t1 = tau(1-u), t2 = tau(1+2u), 4t1+2t2=6tau
% to simulate H = u * Jz \sum_{i<j} (S^y_i S^y_j -1/2 S^x_i S^x_j -1/2
% S^z_i S^z_j) / |i-j|^3 + hz* \sum S^z_j
% 
% general case: u_list = [u v w], (u-w)XX+(v-u)YY+(w-v)ZZ, see MBL PRL
% Chao Yin 6/13
switch nargin
    case 6
        hz = 0;
end
alpha = 3;
cplist = Jz * [1:N_atom-1].^(-alpha) ;
dipolar_z = [-1/2 -1/2 1];

H_nature = OperatorClass(N_atom,dipolar_z,1/4,bc,cplist);
H_nature.symmetrize(sym);

phase = [0 90 90 0 0 90 90 0];
phase = [phase, 180+ phase];
angle = pi/2;
 
u = u_list(1); v = u_list(2); w = u_list(3);
t1 = tau* (1-v+w);
t2 = tau* (1-u+v);
t3 = tau* (1+u-w);
delay = [t1 repmat([t2 2*t3 t2 2*t1], [1,3]) t2 2*t3 t2 t1];
U_real = UFloquet(phase,angle,delay,H_nature);

temp = H2U( OperatorClass(N_atom,'z',1/2) ,hz* 24* tau);
temp.symmetrize(sym);
U_real = temp * U_real;
end

