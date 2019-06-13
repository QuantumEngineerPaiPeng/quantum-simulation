N_atom=3;
tol=1e-4;
maxiter=100;
bc='o';% 'p' for periodic boundary condition; 'o' for open boundary condition
g=0.5;
% cplist=[1,1/8,1/27];
cplist=1;
delta=0.05;
N=1; % number of steps
tlist=(1:N)*delta;

Sigma_i=IndivPauliSparse(N_atom);
Sigma_x=Sigma_i{1}{1};
Sigma_y=Sigma_i{2}{1};
Sigma_z=Sigma_i{3}{1};
for p=2:N_atom
    Sigma_x=Sigma_x+Sigma_i{1}{p};
    Sigma_y=Sigma_y+Sigma_i{2}{p};
    Sigma_z=Sigma_z+Sigma_i{3}{p};
end
o_0=Sigma_z;

H_int=sparse(2^N_atom,2^N_atom);

for pp=1:length(cplist)
    for p=1:N_atom
        if bc=='p'
            q=mod(p+pp-1,N_atom)+1;
        else
            q=p+pp;
            if q>N_atom
                continue
            end
        end
        H_int=H_int+(-1)*cplist(pp)*...
            (Sigma_i{1}{p}*Sigma_i{1}{q}-...
            (Sigma_i{3}{p}*Sigma_i{3}{q}+Sigma_i{2}{p}*Sigma_i{2}{q})/2);
        %         H_int=H_int+(-1)*cplist(pp)*...
        %             (Sigma_i{1}{p}*Sigma_i{1}{q});
        %     H_int=H_int+...
        %         J*(Sigma_i{2}{p}*Sigma_i{2}{p+1}+Sigma_i{3}{p}*Sigma_i{3}{p+1}+Sigma_i{1}{p}*Sigma_i{1}{p+1});
        
    end
end
H=H_int+g*Sigma_z;
tic
[V,D]=eig(full(H));
toc

fuck=mfilename