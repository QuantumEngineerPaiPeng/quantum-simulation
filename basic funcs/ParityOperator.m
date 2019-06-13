% 20180423
% generate parity operator for given number of atoms
function parity=ParityOperator(N_atom)
Sigma_i=IndivPauliSparse(N_atom);
parity=speye(2^N_atom);
for p=1:N_atom
    q=N_atom+1-p;
    if p>=q
        break
    end
    parity=parity*(Sigma_i{1}{p}*Sigma_i{1}{q}+Sigma_i{2}{p}*Sigma_i{2}{q}+...
        Sigma_i{3}{p}*Sigma_i{3}{q}+speye(2^N_atom))/2;
end
end