## Created by Xiaoyang Huang, Apr. 16


import gym
from gym import error, spaces, utils
from gym.utils import seeding
import numpy as np
from scipy.linalg import schur, logm, expm
import operator


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





class XyEnv(gym.Env):
    metadata = {'render.modes': ['human']}

    def __init__(self):
        self.maxTime=None      # total number of pulses in single sequence
        self.nSpin=None
        self.target=None     # target Hamiltonian
        self.unitary=None                      
        self.U0=None          # orginal Hamiltonian
        self.pw=0
        self.H0=None                       # original Hamiltonian
        self.tot_time=0          # total time in one period

        self.state=None
        self.frame=np.eye(2) # Toggling frame unitary

    def setParam(self,maxTime, nSpin,pw):
        self.maxTime=maxTime
        self.action_space=spaces.Box(low=-1, high=1,shape=(10,), dtype=np.float32) ## action 5 for {d,x,y,-x,-y} 5 for corresponding delay or angle shift
                                                      
        self.observation_space = spaces.Box(low=-1, high=1,shape=(25,), dtype=np.float32) ## record ordered action as state with last one the step
        self.nSpin=nSpin
        self.unitary=np.eye(2**nSpin) # Unitary operator
        self.pw=pw

    def setTarget(self,target):
        if self.nSpin!=None:
            sh=target.shape
            if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
                raise ValueError('target must be a 2^nSpin-by-2^nSpin matrix.')
        else:
            raise ValueError('setParam before setTarget')
        self.target=target


    def seed(self, seed=None):
            self.np_random, seed = seeding.np_random(seed)
            return [seed]

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
    def getFidelity(self,n):
        D, V=schur(self.unitary)
        dd=np.diag(np.exp(np.log(np.diag(D))/self.tot_time)) 
    #     dd=np.diag(np.exp(np.log(np.diag(D)))) 

        U=np.dot(np.dot(V,dd),np.transpose(np.conjugate(V))) # unitary operator in unit time
        return np.abs(np.sum(U*np.transpose(np.conjugate(self.target)))/2**self.nSpin)


# # get Tr(self.unitary^n*self.target')
#     def get_n_Fidelity(self,n):
#         D, V=schur(self.unitary)
#     #     dd=np.diag(np.exp(np.log(np.diag(D))/self.state['n'])) 
#         dd=np.diag(np.exp(np.log(np.diag(D)))) 
#         dd_n=dd
#         D1, V1=schur(self.target)
#         dd1=np.diag(np.exp(np.log(np.diag(D1)))) 
#         dd_n1=dd1
#         for pp in range(n):
#             dd_n=np.dot(dd,dd_n)
#             dd_n1=np.dot(dd1,dd_n1)
#         U=np.dot(np.dot(V,dd_n),np.transpose(np.conjugate(V))) # unitary operator in nT time
#         U1=np.dot(np.dot(V1,dd_n1),np.transpose(np.conjugate(V1)))
#         return np.abs(np.sum(U*np.transpose(np.conjugate(U1)))/2**self.nSpin)
    

# get u_frame^{\kron^nSpin}
    def u_frame_nSpin(self, u_frame):
        temp=1
#         print(self.nSpin)
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
    
    
    
    def step(self, action,stepp):
        pw=self.pw

        action = np.clip(action, -1, 1)
        index, value = max(enumerate(action[:5]), key=operator.itemgetter(1)) # get the max action

        phi=action[5+index]*1/2/np.pi+1/2/np.pi+np.pi/2
             
#         phi=action[5+index]*1/np.pi+np.pi/2
        
        self.state[stepp]=index*1.0/4     # record action
        self.state[stepp+12]=action[index+5] # record delay or angle shift
        self.state[24]=stepp/11          # record step


                
        tauu=1
       
        self.tot_time=self.tot_time+tauu+pw*1

        field=self.H_pulse_nSpin(index)
        
        unitary_n=np.dot(expm(-1j*self.H0*tauu),expm(-1j*(self.H0*pw+phi/2*field)))
        self.unitary=np.dot(unitary_n,self.unitary)
        
        if stepp>10:
            reward=-np.log(1-self.getFidelity(stepp+1))
            return self.state, reward, True, {"frame":self.frame}

#         unitary_n=np.dot(np.dot(u_frame,self.U0),np.transpose(np.conjugate(u_frame)))
#         self.unitary=np.dot(unitary_n,self.unitary)

        return self.state, 0.00, False, {"frame":self.frame}      


  
    def reset(self):
        self.tot_time=0
        self.state=np.array([-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.])
        if self.nSpin==None:
            self.unitary=None # Unitary operator
        else:
            self.unitary=np.eye(2**self.nSpin)
        self.frame=np.eye(2)
#         return self.state, {"frame":self.frame}
        return self.state

    def render(self, mode='human', close=False):
        pass
