%20190325
% compute MQC then OTOC in the prethermal paper
% confirm validity of MQC code
% H = Hz + Hx
% Hz = Jz \sum_{i<j} (S^z_i S^z_j -1/2 ... - 1/2 ... ) / |i-j|^alpha + hz \sum S^z_j, 
% Hx = hx \sum S^x_j
% Written using OperatorClass
% Created by Chao Yin
% last change: 3/25

%% runtime parameters
N_atom = 10;

hx_list = [0.25 0.5 1];
Jz = 1;
hz = 0;

hdtau = 0.01;
htau_list = [hdtau];

t_end = 8;
t_intervals = {[1, t_end]}; % negative number triggers ED for infinite long time

load_save = false;
ising = false;

profile on
profile clear

%% use function MQC, C_YZ in paper = C_ZX here
q2_list = [ [0:N_atom], [N_atom-1:-1:1] ].^2;

for hx_index = 1: length(hx_list)
    [I_qtz, I_qtx] = MQC( N_atom, hx_list(hx_index), Jz, hz ...
        ,hdtau./ hx, htau_list./ hx, t_intervals, load_save, ising);
    if hx_index ==1
        C_ZX = sum( I_qtx{1,1} .* q2_list', 1 );
    else
        C_ZX(end+1,:) = sum( I_qtx{1,1} .* q2_list', 1 );
    end
end


%% plot C_ZX
figure(51);
clf;
legends = cell(0);
for hx_index = 1: length(hx_list)
    legends{end+1} = ['$\frac{h_x}{J_z}=',mat2str(hx_list(hx_index) ),'$'];
end
x_data = linspace(t_intervals{1}(1), t_intervals{1}(2) , size( I_qtz{1,1},2) );

ycPlot(x_data, 0.4* C_ZX, '$J_z t$', '$0.4 C_{ZX}$', ...,
    ['dipolar z,$\frac{h_z}{J_z}=',mat2str(hz) ...
    , ' ,N_{atom}=',mat2str(N_atom),' ,h_x \tau=', mat2str(htau_list(1)) ,'$'], legends);




