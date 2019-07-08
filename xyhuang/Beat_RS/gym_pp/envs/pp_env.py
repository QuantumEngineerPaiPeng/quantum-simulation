import gym
from gym import error, spaces, utils
from gym.utils import seeding
import numpy as np
from scipy.linalg import schur, logm, expm



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
            H=H+cplist[pp]*(model[0]*Kron2body(N_atom,1,0,p,q)
                            +model[1]*Kron2body(N_atom,2,0,p,q)
                            +model[2]*Kron2body(N_atom,3,0,p,q))
    if np.max(np.abs(np.imag(H)))<1e-10:                                         #why?
        H=np.real(H)
    return H





class PpEnv(gym.Env):
  metadata = {'render.modes': ['human']}

  def __init__(self):
    self.maxTime=None
    self.nSpin=None
    self.target=None
    self.action_space=spaces.Discrete(6) # delay, x, y, -x, -y, terminate
    self.observation_space=None
    self.state={"pp": 0, "n": 0}
    self.pulses = {"pp": [], "n": 0}
    self.U0=None
    self.pw=None                       # pulse width
    self.H0=None                       # original Hamiltonian
    self.Htarget=None

    self.unitary=None
    self.frame=np.eye(2) # Toggling frame unitary

  def setParam(self, maxTime, nSpin,pw):
    self.maxTime=maxTime
    self.observation_space= spaces.Dict({"pp": spaces.Discrete(5**maxTime), "n": spaces.Discrete(maxTime)})
    self.nSpin=nSpin
    self.unitary=np.eye(2**nSpin) # Unitary operator
    self.pw=pw

  def setTarget(self,target, false_frame=0):
    if self.nSpin!=None:
      sh=target.shape
      if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
        raise ValueError('target must be a 2^nSpin-by-2^nSpin matrix.')
    else:
      raise ValueError('setParam before setTarget')
    self.target=target
    self.false_frame = false_frame
    
    
  def setTargetH(self,target):
        if self.nSpin!=None:
            sh=target.shape
            if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
                raise ValueError('target must be a 2^nSpin-by-2^nSpin matrix.')
        else:
            raise ValueError('setParam before setTarget')
        self.Htarget=target
        
        
        

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
    
  def setH0(self,H0):
        if self.nSpin!=None:
              sh=H0.shape
              if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
                  raise ValueError('H0 must be a 2^nSpin-by-2^nSpin matrix.')
        else:
             raise ValueError('setParam before setH0')
        self.H0=H0

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

  def H_pulse_nSpin(self, index):
        n=index
        if n==0:
            temp=np.eye(2**self.nSpin)
        elif n==1:
            temp=Hamiltonian(self.nSpin,'p',[1],[1,0,0])
        elif n==2:
            temp=Hamiltonian(self.nSpin,'p',[1],[0,1,0])
        elif n==3:
            temp=Hamiltonian(self.nSpin,'p',[1],[-1,0,0])
        elif n==4:
            temp=Hamiltonian(self.nSpin,'p',[1],[0,-1,0])
        return temp
    
  def step(self, action):
    # update self.pulses
    self.pulses['pp'].append(action)
    self.pulses['n'] += 1
    
    if action==5:
      if np.max(np.abs(self.frame-np.eye(2)))<1e-10:
        
        reward=-np.log(1-self.getFidelity())
#         reward=reward_fake
#         if reward<24:
#             reward=self.false_frame
        # reward = 0 if reward<11 else reward
      else:
        reward=self.false_frame
        
      return self.state, reward, True, {"frame":self.frame}
    
    pp=self.state["pp"]+action*5**self.state["n"]
    self.state["pp"]=pp
    self.state["n"]=self.state["n"]+1

    pw=self.pw
    field=self.H_pulse_nSpin(action)
    tauu=1-pw
    phi=np.pi/2
    unitary_n=np.dot(expm(-1j*self.H0*tauu),expm(-1j*(self.H0*pw-phi/2*field)))
    self.unitary=np.dot(unitary_n,self.unitary)
    
    
    if action==0:
      u_pulse=np.eye(2)
    elif action>0 and action<3:
      u_pulse=1/np.sqrt(2)*(np.eye(2)-1j*self.Pauli(action))
    elif action<5:
      u_pulse=1/np.sqrt(2)*(np.eye(2)+1j*self.Pauli(action-2))
    else:
      raise ValueError('Action must be integer from 0 to 5.')
    self.frame=np.dot(self.frame,u_pulse)

    
    
    
    

    return self.state, -0.00, False, {"frame":self.frame}
  
  def reset(self):
    self.state={"pp": 0, "n": 0}
    self.pulses={"pp": [], "n": 0}
    # self.U0=None
    if self.nSpin==None:
      self.unitary=None # Unitary operator
    else:
      self.unitary=np.eye(2**self.nSpin)
    self.frame=np.eye(2)
    return self.state, {"frame":self.frame}
    
  def render(self, mode='human', close=False):
    pass
