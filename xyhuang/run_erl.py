import numpy as np, os, time, sys, random
from core import mod_neuro_evo as utils_ne
from core import mod_utils as utils
import gym, torch
from core import replay_memory
from core import ddpg as ddpg
import argparse
from scipy.linalg import expm, logm
import gym_xy
import math
import matplotlib.pyplot as plt
from copy import deepcopy as dcp
import timeit

render = False
# parser = argparse.ArgumentParser()
# parser.add_argument('-env', help='Environment Choices: (HalfCheetah-v2) (Ant-v2) (Reacher-v2) (Walker2d-v2) (Swimmer-v2) (Hopper-v2)', required=True)
# env_tag = vars(parser.parse_args())['env']

env_tag ='xy-v0'
env = gym.make('xy-v0')


def Pauli(n):
    if n==0:
      return np.eye(2)
    elif n==1:
      return np.array([[0,1],[1,0]])
    elif n==2:
      return np.array([[0,-1j],[1j,0]])
    elif n==3:
      return np.array([[1,0],[0,-1]])
    else:
      raise ValueError('Input must be integer from 0 to 3.')

# returns sigma_a^p*sigma_b^q, with a,b = 1,2,3, p,q being position
def Kron2body(N_atom,a,b,p,q):
    y=1
    for i in range(N_atom):
        if i==p:
            y=np.kron(y,Pauli(a))
        elif i==q:
            y=np.kron(y,Pauli(b))
        else:
            y=np.kron(y,np.eye(2))
    return y

def Hamiltonian(N_atom,bc,cplist,model):
    H=np.zeros((2**N_atom,2**N_atom))
    for pp in range(len(cplist)):
        for p in range(N_atom):
            if bc=='p':
                q=(p+pp+1)%N_atom
            elif bc=='o':
                q=p+pp+1
                if q>=N_atom:
                    continue
            H=H+cplist[pp]*(model[0]*Kron2body(N_atom,1,1,p,q)
                            +model[1]*Kron2body(N_atom,2,2,p,q)
                            +model[2]*Kron2body(N_atom,3,3,p,q))+model[3]*Kron2body(N_atom,3,0,p,q)
    if np.max(np.abs(np.imag(H)))<1e-10:                                         #why?
        H=np.real(H)
    return H


maxTime=12
nSpin=3
pw=1
env.setParam(maxTime,nSpin,pw)


Aim=np.eye(2**env.nSpin)
env.setTarget(Aim)

H=Hamiltonian(nSpin,'p',[1],[-0.5,-0.5,1,0])
J=8.18e-3
t=1
env.setU0(expm(-1j*J*H*t))   
env.setH0(J*H)  



class Parameters:
    def __init__(self):

        #Number of Frames to Run
        self.num_frames = 800000

        #USE CUDA
        self.is_cuda = False; self.is_memory_cuda = False

        #Sunchronization Period
        self.synch_period = 10

        #DDPG params
        self.use_ln = True
        self.gamma = 0.99; self.tau = 0.001
        self.seed = 7
        self.batch_size = 128
        self.buffer_size = 1000000
        self.frac_frames_train = 1.0
        self.use_done_mask = True

        ###### NeuroEvolution Params ########
        #Num of trials
        self.num_evals = 1

        #Elitism Rate
        self.elite_fraction = 0.2


        self.pop_size = 10
        self.crossover_prob = 0.1
        self.mutation_prob = 0.9

        #Save Results
        self.state_dim = None; self.action_dim = None #Simply instantiate them here, will be initialized later
        self.save_foldername = 'R_ERL/'
        if not os.path.exists(self.save_foldername): os.makedirs(self.save_foldername)
            
            
class Agent:
    def __init__(self, args, env):
        self.args = args; self.env = env
        self.evolver = utils_ne.SSNE(self.args)
        self.best_r=0
        self.best_state=[]

        #Init population
        self.pop = []
        for _ in range(args.pop_size):
            self.pop.append(ddpg.Actor(args))

        #Turn off gradients and put in eval mode
        for actor in self.pop: actor.eval()

        #Init RL Agent
        self.rl_agent = ddpg.DDPG(args)
        self.replay_buffer = replay_memory.ReplayMemory(args.buffer_size)
        self.ounoise = ddpg.OUNoise(args.action_dim)

        #Trackers
        self.num_games = 0; self.num_frames = 0; self.gen_frames = None

    def add_experience(self, state, action, next_state, reward, done):
        reward = utils.to_tensor(np.array([reward])).unsqueeze(0)
        if self.args.is_cuda: reward = reward.cuda()
        if self.args.use_done_mask:
            done = utils.to_tensor(np.array([done]).astype('uint8')).unsqueeze(0)
            if self.args.is_cuda: done = done.cuda()
        action = utils.to_tensor(action)
        if self.args.is_cuda: action = action.cuda()
        self.replay_buffer.push(state, action, next_state, reward, done)

    def evaluate(self, net, is_render, is_action_noise=False, store_transition=True):
        total_reward = 0.0

        state = self.env.reset()
        state = utils.to_tensor(state).unsqueeze(0)
#         state = utils.to_tensor(state)
        if self.args.is_cuda: state = state.cuda()
        done = False
        pp=0

        while not done:
            if store_transition: self.num_frames += 1; self.gen_frames += 1
            if render and is_render: self.env.render()
            action = net.forward(state)
            action.clamp(-1,1)
            action = utils.to_numpy(action.cpu())
            if is_action_noise: action += self.ounoise.noise()

            next_state, reward, done, info = self.env.step(action.flatten(),pp)  #Simulate one step in environment
            next_state = utils.to_tensor(next_state).unsqueeze(0)
            if self.args.is_cuda:
                next_state = next_state.cuda()
            total_reward += reward

            if store_transition: self.add_experience(state, action, next_state, reward, done)
            state = next_state
            pp=pp+1
        if store_transition: self.num_games += 1
        if total_reward>self.best_r:
            self.best_r=total_reward
            self.best_state.append(state)
            print(state)

        return total_reward

    def rl_to_evo(self, rl_net, evo_net):
        for target_param, param in zip(evo_net.parameters(), rl_net.parameters()):
            target_param.data.copy_(param.data)

    def train(self):
        self.gen_frames = 0

        ####################### EVOLUTION #####################
        all_fitness = []
        #Evaluate genomes/individuals
        for net in self.pop:
            fitness = 0.0
            for eval in range(self.args.num_evals): fitness += self.evaluate(net, is_render=False, is_action_noise=False)
            all_fitness.append(fitness/self.args.num_evals)

        best_train_fitness = max(all_fitness)
        worst_index = all_fitness.index(min(all_fitness))

        #Validation test
        champ_index = all_fitness.index(max(all_fitness))
        test_score = 0.0
        for eval in range(5): test_score += self.evaluate(self.pop[champ_index], is_render=True, is_action_noise=False, store_transition=False)/5.0

        #NeuroEvolution's probabilistic selection and recombination step
        elite_index = self.evolver.epoch(self.pop, all_fitness)


        ####################### DDPG #########################
        #DDPG Experience Collection
        self.evaluate(self.rl_agent.actor, is_render=False, is_action_noise=True) #Train

        #DDPG learning step
        if len(self.replay_buffer) > self.args.batch_size * 5:
            for _ in range(int(self.gen_frames*self.args.frac_frames_train)):
                transitions = self.replay_buffer.sample(self.args.batch_size)
                batch = replay_memory.Transition(*zip(*transitions))
                self.rl_agent.update_parameters(batch)

            #Synch RL Agent to NE
            if self.num_games % self.args.synch_period == 0:
                self.rl_to_evo(self.rl_agent.actor, self.pop[worst_index])
                self.evolver.rl_policy = worst_index
                print('Synch from RL --> Nevo')

        return best_train_fitness, test_score, elite_index
    
    
start = timeit.default_timer()    
average_score=[]
frame=[]
best_fid=0
best_state=None
parameters = Parameters()  # Create the Parameters class
tracker = utils.Tracker(parameters, ['erl'], '_score.csv')  # Initiate tracker
frame_tracker = utils.Tracker(parameters, ['frame_erl'], '_score.csv')  # Initiate tracker
time_tracker = utils.Tracker(parameters, ['time_erl'], '_score.csv')

#Create Env
#     env = utils.NormalizedActions(gym.make(env_tag))
parameters.action_dim = env.action_space.shape[0]
parameters.state_dim = env.observation_space.shape[0]

    #Seed
env.seed(parameters.seed);
torch.manual_seed(parameters.seed); np.random.seed(parameters.seed); random.seed(parameters.seed)

    #Create Agent
agent = Agent(parameters, env)
print('Running', env_tag, ' State_dim:', parameters.state_dim, ' Action_dim:', parameters.action_dim)

next_save = 100; time_start = time.time()

##### loading checkpoint
agent.pop[0].load_state_dict(torch.load(parameters.save_foldername + 'evo_net'))
agent.pop[0].eval()
while agent.num_frames <= parameters.num_frames:
    best_train_fitness, erl_score, elite_index = agent.train()
    if best_train_fitness>best_fid:
        best_fid=best_train_fitness
    average_score.append(tracker.all_tracker[0][1])
    frame.append(agent.num_frames)
            
    print('#Games:', agent.num_games, '#Frames:', agent.num_frames, ' Epoch_Max:', '%.2f'%best_train_fitness if best_train_fitness != None else None, ' Test_Score:','%.2f'%erl_score if erl_score != None else None, ' Avg:','%.2f'%tracker.all_tracker[0][1], 'ENV '+env_tag)
    print('RL Selection Rate: Elite/Selected/Discarded', '%.2f'%(agent.evolver.selection_stats['elite']/agent.evolver.selection_stats['total']),
                                                             '%.2f' % (agent.evolver.selection_stats['selected'] / agent.evolver.selection_stats['total']),
                                                              '%.2f' % (agent.evolver.selection_stats['discarded'] / agent.evolver.selection_stats['total']))
    print()
    tracker.update([erl_score], agent.num_games)
    frame_tracker.update([erl_score], agent.num_frames)
    time_tracker.update([erl_score], time.time()-time_start)

    #Save Policy
    if agent.num_games > next_save:
        next_save += 100
        if elite_index != None: torch.save(agent.pop[elite_index].state_dict(), parameters.save_foldername + 'evo_net')
        print("Progress Saved")
            
stop = timeit.default_timer()

print('Time: ', stop - start)          
            

print('Parameters.num_frames: ',parameters.num_frames)
print('Parameters.synch_period: ',parameters.synch_period)
print('Parameters.gamma: ',parameters.gamma)
print('Parameters.tau: ',parameters.tau)
print('Parameters.seed: ',parameters.seed)
print('Parameters.batch_size: ',parameters.batch_size)
print('Parameters.buffer_size: ',parameters.buffer_size)
print('Parameters.frac_frames_train: ',parameters.frac_frames_train)
print('Parameters.num_evals: ',parameters.num_evals)
print('Parameters.pop_size: ',parameters.pop_size)
print('Parameters.crossover_prob: ',parameters.crossover_prob)
print('Parameters.mutation_prob: ',parameters.mutation_prob)



average_score[0]=average_score[1]

np.savetxt('ES_DDPG_pw_12_angleshift_gradinv_period10_gamma0.99_taue-3_seed7_batchsize128_buffersizee6_fracframestrain1_numevals1_elitefraction0.2_popsize10_crossoverprob0.1_mutationprob0.9.csv',np.array([frame,  average_score]).T,fmt=['%.7f','%.7f'],delimiter=',',header="co,rew") 

fig, ax = plt.subplots(figsize=(5,4),tight_layout=True)

ax.plot(frame,  average_score)


# ax.set_xscale('log')

ax.set(xlabel='Frame', ylabel='avrReward',
       title='Jt=8.18e-3')        ### The 12 pulse sequence is compared with the 6 pulse sequence
                                   ### We regard the 6 pulse sequence is T long, the 12 pulse sequence is 2T long 

# plt.legend(['Fid_ML12', 'Fid_ML_sym12','Fid_ML6', 'Fid_WHH'], loc='best')

plt.savefig('ES_DDPG_pw_12_angleshift_gradinv_period10_gamma0.99_taue-3_seed7_batchsize128_buffersizee6_fracframestrain1_numevals1_elitefraction0.2_popsize10_crossoverprob0.1_mutationprob0.9.eps', dpi=fig.dpi, bbox_inches='tight')

temp=agent.best_state[len(agent.best_state)-1][0].numpy()
temp=np.append(temp,agent.best_r)

np.savetxt('ES_DDPG_pw_12_angleshift_gradinv.csv',np.array(temp).T,fmt=['%.7f'],delimiter=',',header="rew") 


print(best_fid)
















