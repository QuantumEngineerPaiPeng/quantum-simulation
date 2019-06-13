function S_mt = MQC_S( ...,
    U_t, O0, m_list, t_list, PZ ...,
)
% 2019/3/23
% S_mt = S_m(t) = 2^(-L) tr( O0 U(m,t)' O0 U(m,t) )
% U(m, t) = U_t^(-t) exp(im pi PZ/L) U_t^t
% use OperatorClass
% Chao Yin
% last change 3/23
S_mt = zeros(length(m_list), length(t_list));
if nargin <= 4
else
    for m_index = 1: length(m_list)
        e_P = H2U(PZ, -pi/U_t.L * m_list(m_index));
        parfor t_index = 1: length(t_list)
            U_tt = U_t^t_list(t_index);
            U_mt = U_tt' * e_P * U_tt;
            
            assignin('base','U_M',U_tt'* O0 * U_tt);
            S_mt(m_index, t_index) = trace( O0* U_mt' * O0 * U_mt );
        end
    end
end
S_mt = real( S_mt)/ 2^U_t.L;

end

