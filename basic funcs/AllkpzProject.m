% 20180423
% generate all the pt projection operators for given # of atoms
function y=AllkpzProject(N_atom)

tempy=cell((floor(N_atom/2)+1),1);
parfor k=1:floor(N_atom/2)+1
    temp={kpzProject(N_atom,k-1,-1,-1),kpzProject(N_atom,k-1,1,-1),...
        kpzProject(N_atom,k-1,-1,1),kpzProject(N_atom,k-1,1,1)};
    tempy{k}=temp;
end
y={};
for k=1:floor(N_atom/2)+1
    y{end+1}=tempy{k}{1};
    y{end+1}=tempy{k}{2};
end
for k=1:floor(N_atom/2)+1
    y{end+1}=tempy{k}{3};
    y{end+1}=tempy{k}{4};
end
end