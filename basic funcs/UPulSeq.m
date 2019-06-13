% 20180622
% Generate the unitary oeprators for each pulse in a pulse sequence
% Input: 
% L: system size
% phase: array, giving the phases of pulse sequence in unit of degree
% angle: array or scalar (same for all pulses), giving the angles of pulse
% sequence in unit of rad
% sym: optional, symmtry class
% Output: cell, containing the unitaries
% created by Pai Peng
function y=UPulSeq(L,phase,angle,sym)
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

phase=phase/180*pi;

X=OperatorClass(L,'x',1);
Y=OperatorClass(L,'y',1);
if nargin==4
    X.symmetrize(sym)
    Y.symmetrize(sym)
end

y=cell(1,length(phase));
for p=1:length(phase)
    y{p}=H2U(X*cos(phase(p))+Y*sin(phase(p)),angle(p)/2);
end