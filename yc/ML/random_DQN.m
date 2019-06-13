%% read data from npz
data = unzip('E:\career\junior2\code\RL-qtm-engineering\DQN_T12.npz');
n_run = readNPY('n_run.npy');
best_epi = readNPY('best_episode.npy');

%% histogram
figure(1)
clf;

subplot(1,2,1)
hold on
his_edge = [0 50 100 150 200 250 300];
histogram(best_epi(1,:))%, his_edge)
histogram(best_epi(2,:))%, his_edge)
legend('random','DQN')
title('nEpisode=3e3, lr=0.01, n\_run=100')

subplot(1,2,2)
hold on
his_edge = [0 50 100 150 200 250 300];
histogram(best_epi(1,:), his_edge)
histogram(best_epi(2,:), his_edge)
legend('random','DQN')
title('nEpisode=3e3, lr=0.01, n\_run=100')

%% deleta npy
delete('n_run.npy','best_episode.npy');