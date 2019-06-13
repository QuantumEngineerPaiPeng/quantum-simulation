%20190327
% find out whats going on at hx* tau = pi/2 using map to 2d ising
% Created by Chao Yin
% last change: 3/27

%% runtime parameters
N = 1;

hx = 1;
Jz = 1;
hz = 1;

htau = pi/2;
tau = htau/ hx;

t_end = 4;
M = t_end-1;

profile on
profile clear

%% O(0) = Z
Z = zeros( 2^N, 2^N );
for r_index = 1: 2^N
    ss = dec2bin(r_index-1, N)* 2 - 97;
    Z(r_index, r_index) = sum(ss);
end % for r_index = 1: 2^N

%% U_t O(0) U_-t part in S_m
K1 = i* tau* Jz;
K2 = i*pi/4;
T1 = Ising2d_T( N, K1/4, K2/4 );
T2 = Ising2d_T( N, -K1/4, K2'/4 );

UOU = T1^(M+1) * Z * T2^(M+1);
abs_UOU = abs(UOU);

figure(55)
pcolor( enlarge_pcolor( abs_UOU) )
colorbar()

%% calculate S_m
S_m = zeros(2*N, 1);
for m = 0: 2*N-1
    S_m(m+1) = trace(Z* T1^(M+1)* expm(i*m*pi/N*Z/2) * T2^(M+1)* Z* T1^(M+1)* expm(-i*m*pi/N*Z/2) * T2^(M+1))/ 2^N;
end

%% sum over odd boundarys and even boundaries
sum_ = [0 0];
[zero_n, max_n] = deal(0);
abs_UOU = abs( U_M.matrix{1} );
N = N_atom;  
for r_index = 1: 2^N
    ss = dec2bin(r_index-1, N)-48; % 1   or 0
    for c_index = 1: 2^N
        sp = dec2bin(c_index-1, N)-48;
        if abs_UOU(r_index, c_index) < 1e-10
            zero_n = zero_n+1;
        elseif abs(0.25*sqrt(2)- abs_UOU(r_index, c_index)) < 1e-10
            max_n = max_n+1;
        end
        
        if mod( sum(ss)+sum(sp), 2) == 0
            %[r_index, c_index]
            sum_(1) = sum_(1) + abs_UOU(r_index, c_index);
        elseif mod( sum(ss)+sum(sp), 2) == 1
            sum_(2) = sum_(2) + abs_UOU(r_index, c_index);
        else
            error('sum(ss)+sum(sp) \neq int')
        end
    end
end % for r_index = 1: 2^N

sum_
zero_n


