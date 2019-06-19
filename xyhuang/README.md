# Hamiltonian engineering by RL
Code referring to  "Evolution-Guided Policy Gradients in Reinforcement Learning" accepted at NIPS 2018\
Gradient inverse referring to "DEEP REINFORCEMENT LEARNING IN PARAMETERIZED ACTION SPACE"\

###### Dependencies #######
Python 3.5.6 \
Pytorch 0.3.1.post3 \
Numpy 1.15.2 \
Fastrand from https://github.com/lemire/fastrand \
Gym 0.10.5 \


#### Hyperparameters #### 
self.num_frames = 5000000\
self.is_cuda = False; self.is_memory_cuda = False\
self.synch_period = 10\
self.use_ln = True\
self.gamma = 0.99; self.tau = 0.001\
self.seed = 7\
self.batch_size = 128\
self.buffer_size = 1000000\
self.frac_frames_train = 1.0\
self.use_done_mask = True\

self.num_evals = 1

        
self.elite_fraction = 0.2


self.pop_size = 10\
self.crossover_prob = 0.1\
self.mutation_prob = 0.9\

actor learning rate: lr=0.5e-5\
critic learning rate: lr=0.5e-4\
for 6 pulses neurons are 128 for actor 300 for critic\
for 12 pulses neurons are 256 for actor 400 for critic\
#### To run on cluster #### 
loading checkpoint first on your laptop

sbatch run.sh


