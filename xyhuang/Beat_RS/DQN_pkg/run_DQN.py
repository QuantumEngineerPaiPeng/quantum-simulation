#!/usr/bin/env python
# coding: utf-8

# # Reinforcement learning for Hamiltonian engineering

# ## Set environment

# In[19]:


import os, subprocess, time, signal
import gym
from gym import error, spaces
from gym import utils
from gym.utils import seeding
import numpy as np
from scipy.linalg import expm
from scipy import sparse
import gym_pp

from DQN import DeepQNetwork

# 10-base to arbitrary base
def dec2base(num,base,length):
    s=''
    if num>base**length-1:
        raise ValueError('Input number exceeds the maximum number allowed by the length')
    for i in range(length):
        s=s+chr(ord('0')+int(num/(base**(length-1-i))))
        num=num-int(num/(base**(length-1-i)))*(base**(length-1-i))
    return s

# given a state dict, return a str that specifies all pulses
def getPPstr(state,actionDict=['d','x','y','-x','-y']):
    num=state['pp']
    base=len(actionDict)
    length=state['n']
    rawPulseStr=dec2base(num,base,length)
    pulseStr=''
    for p in rawPulseStr:
        pulseStr=actionDict[int(p)]+','+pulseStr
    pulseStr=pulseStr[0:-1]
    return pulseStr

env = gym.make('pp-v0')

observation = env.reset()


# ## Set parameters

# In[20]:


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
                            +model[2]*Kron2body(N_atom,3,3,p,q))
    if np.max(np.abs(np.imag(H)))<1e-10:
        H=np.real(H)
    return H


# In[21]:


maxTime=8
nSpin=3
env.setParam(maxTime,nSpin)

env.setTarget(np.eye(2**env.nSpin))
H=Hamiltonian(nSpin,'p',[1],[-0.5,-0.5,1])
J=8.18e-3
t=1
env.setU0(expm(-1j*J*H*t))


# ## Q learning

# In[22]:


nEpisodes=100
alpha=0.1 # gradient step
gamma=0.9 # discount factor
# beta=0.2*np.ones(nEpisodes) # inverse temperature
#beta=np.logspace(-3.0, -.5, num=nEpisodes) # inverse temperature

#qTable=np.zeros((5**maxTime,6))
# qTable=sparse.csr_matrix((5**maxTime,6))

# how many actions are available (maxTime is not considered)
def getAvailableAction(state,frame):
    if np.max(np.abs(frame-np.eye(2)))<1e-10 and state['n']!=0:
        return 6
    else:
        return 5

# get the next action
def getAction(state,frame,beta):
    if state['n']==maxTime:
        return 5
    nAction=getAvailableAction(state,frame)
    prob=np.exp(beta*qTable[state['pp'],0:nAction])
    prob=prob/np.sum(prob)
    action=np.random.choice(nAction,1,p=prob)
    return action[0]

# DQN
def action_transform(pp_obser):
    return [pp_obser['pp'], pp_obser['n']]

env = env.unwrapped

print(env.action_space)
print(env.observation_space)
#print(env.observation_space.high)
#print(env.observation_space.low)

RL = DeepQNetwork(n_actions=env.action_space.n,
                  n_features=2,
                  learning_rate=alpha, reward_decay=gamma, e_greedy=0.9,
                  replace_target_iter=100, memory_size=2000,
                  e_greedy_increment=0.01,)
total_steps = 0




best_reward=0
best_pp=None
reward_list=[]
for episode in range(nEpisodes):
    # DQN
    observation, info = env.reset()
    # frame=info["frame"]
    ep_r = 0
    while True:
        env.render()

        action = RL.choose_action(action_transform(observation))

        observation_, reward, done, info = env.step(action)

        
        RL.store_transition(action_transform(observation), action, reward, action_transform(observation_))

        ep_r += reward
        if total_steps > 1000:
            RL.learn()

        if done:
            print('episode: ', episode,
                  'ep_r: ', round(ep_r, 2),
                  ' epsilon: ', round(RL.epsilon, 2))
            break

        observation = observation_
        total_steps += 1



    if reward>best_reward:
        best_reward=reward
        best_pp=getPPstr(env.state)
        print('Episode: '+str(episode))
        print(getPPstr(env.state))
        print('Fidelity: '+str(env.getFidelity())+', Reward: '+str(reward))
    reward_list.append(reward)
    if (episode+1)%1000==0 and episode!=0:
        recent_reward=reward_list[episode-99:episode+1]
        print('Recent 1000 average reward: '+str(np.mean(recent_reward)))
    

       
print('Best sequece: '+best_pp)


# ## Sequence tests

# In[15]:


H=np.array([[3,0,0,0,0,0,0,0],
     [0,-1,-1,0,-1,0,0,0],
     [0,-1,-1,0,-1,0,0,0],
     [0,0,0,-1,0,-1,-1,0],
     [0,-1,-1,0,-1,0,0,0],
     [0,0,0,-1,0,-1,-1,0],
     [0,0,0,-1,0,-1,-1,0],
     [0,0,0,0,0,0,0,3]])
J=8.18e-3
t=1
env.setU0(expm(-1j*J*H*t*4))


# In[16]:


env.reset()
observation, reward, done, info = env.step(3)

observation, reward, done, info = env.step(4)

observation, reward, done, info = env.step(2)

observation, reward, done, info = env.step(1)

observation, reward, done, info = env.step(2)

observation, reward, done, info = env.step(4)

observation, reward, done, info = env.step(5)
print(getPPstr(env.state))
print(env.frame)
print('Fidelity: '+str(env.getFidelity())+', Reward: '+str(reward))


# In[18]:


# WAHUHA
env.reset()
observation, reward, done, info = env.step(0)

observation, reward, done, info = env.step(1)

observation, reward, done, info = env.step(4)

observation, reward, done, info = env.step(0)

observation, reward, done, info = env.step(2)

observation, reward, done, info = env.step(3)

observation, reward, done, info = env.step(5)

print('WAHUHA: '+getPPstr(env.state))
print(env.frame)
print('Fidelity: '+str(env.getFidelity())+', Reward: '+str(reward))


# In[ ]:





# In[ ]:





# In[ ]:




