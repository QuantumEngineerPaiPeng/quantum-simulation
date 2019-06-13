% 20180911
% function to generate Hamiltonian in free Fermion picture
% input:
% model: 2-element array, giving the coefficient in front of XX and YY
% N_spin: system size
% bc: 'o' or 'p', boundary condition
% Created by Pai Peng
function y=FermionHamiltonian(model,N_spin,bc)
H11=diag(ones(1,N_spin-1),1);
if bc=='p'
    H11(1,N_spin)=-1;
end
H12=H11-transpose(H11);
if bc=='p'
    H12(1,N_spin)=1;
    H12(N_spin,1)=-1;
end
H11=H11+transpose(H11);
XX=-1/2*[H11,H12;
    -H12,-H11];
YY=-1/2*[H11,-H12;
    +H12,-H11];
y=model(1)*XX+model(2)*YY;