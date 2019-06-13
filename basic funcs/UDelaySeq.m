% 20180623
% Generate the unitary oeprators for each delay in a pulse sequence
% Input: 
% H: Hamiltonian, Operator class
% delay: array, giving the delay time
% Output: cell, containing the unitaries
% created by Pai Peng
function y=UDelaySeq(H,delay)
if ~isvector(delay)
    error('Delay must be given in a vector')
end

y=cell(1,length(delay));
for p=1:length(delay)
    y{p}=H2U(H,delay(p));
end