function T = Ising2d_T( ...,
    N, K1, K2 ...,
)
% 2019/3/27
% transfer matrix of 2d Ising, T: 2^N * 2^N
% T(s1...sN, s1' ... sN') = exp( -1/2 H(s1...sN) -1/2 H(s1' ... sN') )
% * exp(-K2* (s1s1'+ ... + sNsN'))
% H(s1...sN) = K1* (s1s2+ ... +sNs1)
% s1 = \pm 1
% within each line: N spins, K1, periodic
% inter line: K2
% Chao Yin
% last change 3/27
T = zeros( 2^N, 2^N );
for r_index = 1: 2^N
    ss = dec2bin(r_index-1, N)* 2 - 97;
    Hss = sum( K1/2 * ss .* circshift(ss,1) );
    for c_index = 1: 2^N
        sp = dec2bin(c_index-1, N)* 2 - 97;
        
        Hsp = sum( K1/2 * sp .* circshift(sp,1) );
        HK2 = sum( K2 * sp .* ss );
        T(r_index, c_index) = exp( -Hss -Hsp -HK2 );
    end
end % for r_index = 1: 2^N

end

