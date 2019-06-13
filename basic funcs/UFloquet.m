% 20180623
% construct Floquet unitary for a given pulse sequence
% Assume the sequence starts with a delay
% Input:
% phase: array, giving the phases of pulse sequence in unit of degree
% angle: array or scalar (same for all pulses), giving the angles of pulse
% sequence in unit of rad
% delay: array, giving the delay time
% H: Hamiltonian, Operator class
% W: optional, array or scalar (same for all pulses), giving the widths of pulse
% No W given assumes delta pulses
% Output: unitary, given as OperatorClass
% created by Pai Peng
function y=UFloquet(phase,angle,delay,H,W)
% phase=[0,90];
% angle=pi;
% delay=[1,1];
% H=H_int;
if ~(length(phase)==length(delay) || (length(phase)+1==length(delay)))
    error('Number of pulses must be equal to or 1 less than the number of delays')
end

switch nargin
    case 4
        Upulse=UPulSeq(H.L,phase,angle,H.sym);
    case 5
        Upulse=UPulSeqFW(phase,angle,W,H);
    otherwise
        error('Invalid number of inputs')
end

Udelay=UDelaySeq(H,delay);

y=Udelay{1}';
for p=1:length(phase)
    y=y*(Upulse{p}');
    if (1+p)<=length(delay)
        y=y*(Udelay{1+p})';
    end
end
y=y';