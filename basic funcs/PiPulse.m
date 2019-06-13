% 20180423
% generate the unitary of pi pulse along given direction
function y=PiPulse(N_atom,direction)
sigmalist={sparse(sigma_x),sparse(sigma_y),sparse(sigma_z)};
y=speye(1);
for p=1:N_atom
    y=kron(y,sigmalist{direction});
end
end