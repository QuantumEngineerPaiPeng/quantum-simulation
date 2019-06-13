% 20190607
% returns sigma_a^p*sigma_b^q*sigma_c^r, with a,b,c = 1,2,3, p,q,r being position
function y=Kron3body(N_atom,a,b,c, p,q,r)
sigmalist={sparse(sigma_x),sparse(sigma_y),sparse(sigma_z)};
y=speye(1);
for i=1:N_atom
    switch i
        case p
            y=kron(y,sigmalist{a});
        case q
            y=kron(y,sigmalist{b});
        case r
            y=kron(y,sigmalist{c});
        otherwise
            y=kron(y,speye(2));
    end
end