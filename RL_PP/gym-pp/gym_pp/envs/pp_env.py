import gym
from gym import error, spaces, utils
from gym.utils import seeding
import numpy as np
from scipy.linalg import schur

class PpEnv(gym.Env):
  metadata = {'render.modes': ['human']}

  def __init__(self):
    self.maxTime=None
    self.nSpin=None
    self.target=None
    self.action_space=spaces.Discrete(6) # delay, x, y, -x, -y, terminate
    self.observation_space=None
    self.state={"pp": 0, "n": 0}
    self.U0=None

    self.unitary=None
    self.frame=np.eye(2) # Toggling frame unitary

  def setParam(self, maxTime, nSpin):
    self.maxTime=maxTime
    self.observation_space= spaces.Dict({"pp": spaces.Discrete(5**maxTime), "n": spaces.Discrete(maxTime)})
    self.nSpin=nSpin
    self.unitary=np.eye(2**nSpin) # Unitary operator

  def setTarget(self,target):
    if self.nSpin!=None:
      sh=target.shape
      if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
        raise ValueError('target must be a 2^nSpin-by-2^nSpin matrix.')
    else:
      raise ValueError('setParam before setTarget')
    self.target=target

  def Pauli(self, n):
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

  def setU0(self,U0):
    if self.nSpin!=None:
      sh=U0.shape
      if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
        raise ValueError('U0 must be a 2^nSpin-by-2^nSpin matrix.')
    else:
      raise ValueError('setParam before setU0')
    self.U0=U0

# get Tr(self.unitary*self.target')
  def getFidelity(self):
    D, V=schur(self.unitary)
    dd=np.diag(np.exp(np.log(np.diag(D))/self.state['n'])) 
    U=np.dot(np.dot(V,dd),np.transpose(np.conjugate(V))) # unitary operator in unit time
    return np.abs(np.sum(U*np.transpose(np.conjugate(self.target)))/2**self.nSpin)
    

# get u_frame^{\kron^nSpin}
  def u_frame_nSpin(self, u_frame):
    temp=1
    for p in range(self.nSpin):
      temp=np.kron(temp,u_frame)
    return temp
    
  def step(self, action):
    if action==5:
      if np.max(np.abs(self.frame-np.eye(2)))<1e-10:
        reward=-np.log(1-self.getFidelity())
      else:
        reward=-1000
      return self.state, reward, True, {"frame":self.frame}
    
    pp=self.state["pp"]+action*5**self.state["n"]
    self.state["pp"]=pp
    self.state["n"]=self.state["n"]+1
    if action==0:
      u_pulse=np.eye(2)
    elif action>0 and action<3:
      u_pulse=1/np.sqrt(2)*(np.eye(2)-1j*self.Pauli(action))
    elif action<5:
      u_pulse=1/np.sqrt(2)*(np.eye(2)+1j*self.Pauli(action-2))
    else:
      raise ValueError('Action must be integer from 0 to 5.')
    self.frame=np.dot(self.frame,u_pulse)
    u_frame=self.u_frame_nSpin(self.frame)
    unitary_n=np.dot(np.dot(u_frame,self.U0),np.transpose(np.conjugate(u_frame)))
    self.unitary=np.dot(unitary_n,self.unitary)
    return self.state, -0.01, False, {"frame":self.frame}
  
  def reset(self):
    self.state={"pp": 0, "n": 0}
    # self.U0=None
    if self.nSpin==None:
      self.unitary=None # Unitary operator
    else:
      self.unitary=np.eye(2**self.nSpin)
    self.frame=np.eye(2)
    return self.state, {"frame":self.frame}
    
  def render(self, mode='human', close=False):
    pass
