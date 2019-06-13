% 20180423
% returns sigma_a^p*sigma_b^q, with a,b = 1,2,3, p,q being position
function y=Kron2body(N_atom,a,b,p,q)
sigmalist={sparse(sigma_x),sparse(sigma_y),sparse(sigma_z)};
y=speye(1);
for i=1:N_atom
    switch i
        case p
            y=kron(y,sigmalist{a});
        case q
            y=kron(y,sigmalist{b});
        otherwise
            y=kron(y,speye(2));
    end
end