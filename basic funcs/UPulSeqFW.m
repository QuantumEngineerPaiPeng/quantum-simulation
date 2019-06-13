% 20180622
% Generate the unitary oeprators for each pulse in a pulse sequence,
% considering the finite width of the pulses.
% Input: 
% L: system size
% phase: array, giving the phases of pulse sequence in unit of degree
% angle: array or scalar (same for all pulses), giving the angles of pulse
% sequence in unit of rad
% W: array or scalar (same for all pulses), giving the widths of pulse
% H: Hamiltonian, OperatorClass
% Output: cell, containing the unitaries
% created by Pai Peng
function y=UPulSeqFW(phase,angle,W,H)
L=H.L;
% phase=[0,0];
% angle=pi;
% W=0.001;
% H=H_int;
if ~isvector(phase)
    error('Phase must be given in a vector')
end
if isscalar(angle)
    angle=angle*ones(size(phase));
else
    if length(phase)~=length(angle)
        error('Angle must be a scalar or a vector of the same length as phase')
    end
end

if isscalar(W)
    W=W*ones(size(phase));
else
    if length(W)~=length(phase)
        error('Width must be a scalar or a vector of the same length as phase')
    end
end

phase=phase/180*pi;

X=OperatorClass(L,'x',1);
Y=OperatorClass(L,'y',1);

if ~isempty(H.sym)
    X.symmetrize(H.sym)
    Y.symmetrize(H.sym)
end

y=cell(1,length(phase));
for p=1:length(phase)
    y{p}=H2U((X*cos(phase(p))+Y*sin(phase(p)))*(angle(p)/2/W(p))+H,W(p));
end