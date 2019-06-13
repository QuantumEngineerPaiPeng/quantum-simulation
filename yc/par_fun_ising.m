function Z = par_fun_ising( ...,
    N1, N2, K1, K2, h, method ...,
)
% 2019/4/1
% compute partition function Z of 2d ising model
% Z = \sum e^{-H}
% H = K1 \sigma_j,k \sigma_j+1,k + K2 \sigma_j,k \sigma_j,k+1 + h
% \sigma_j,k
% in our case, K_1 = -i\tau Jz, K_2= i\pi/4 - 1/2 ln(tan(hx \tau)), 
% h = -i\tau hz
% Created by Chao Yin
% last change: 4/1
if nargin == 5 | strcmp( method{1}, 'enumerate')
    Z = 0;
    for state_i = 0: 2^(N1* N2)-1
        state = reshape( dec2bin(state_i, N1*N2)*2 -97 , N1, N2 );
        energy = 0;
        for row_i = 1: N1
            for col_i = 1: N2
                energy = energy + state(row_i, col_i)* ( K1* state( mod(row_i, N1)+1 , col_i) ...
                    + K2*  state( row_i, mod(col_i, N2)+1) + h);
            end
        end
        Z = Z + exp(-energy);
    end % for state_i = 0: 2^(N1* N2)-1
    Z = Z / 2^(N1* N2);
elseif strcmp( method{1}, 'monte')
    Z = 0;
    for state_i = 1: method{2}
        state = 2*randi(2, N1, N2)-3;
        energy = 0;
        for row_i = 1: N1
            for col_i = 1: N2
                energy = energy + state(row_i, col_i)* ( K1* state( mod(row_i, N1)+1 , col_i) ...
                    + K2*  state( row_i, mod(col_i, N2)+1) + h);
            end
        end
        Z = Z + exp(-energy);
    end
    Z = Z/ method{2};
end % if nargin == 5

end

