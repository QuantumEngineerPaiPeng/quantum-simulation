% 20180911
% function to generate disorder Hamiltonian in free Fermion picture
% input:
% disorder: disorder realization
% Created by Pai Peng
function y=FermionDisorderHam(disorder)
L=length(disorder);
D=diag(disorder);
y=[D,zeros(L);zeros(L),-D];