## Created by Xiaoyang Huang, Apr. 16
# We continuously apply pulses , and try to maintain high Fidelity at the end of each sequence
# The Bang-Bang pi/2 pulses are in continuous action space with continuous delay time [min_interv, max_interv]. Using 1d real axis to distinguish different types of pulses
# 

import gym
from gym import error, spaces, utils
from gym.utils import seeding
import numpy as np
from scipy.linalg import schur, logm, expm
import operator
import scipy.integrate as integrate

def commute(a,b):
    return np.dot(a,b)-np.dot(b,a)

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
        self.targetH=None     # target Hamiltonian
        self.unitary=None                      
        self.U0=None          # orginal Hamiltonian
        self.min_delay=None
        self.max_delay=None
        self.pw=None                       # pulse width
        self.H0=None                       # original Hamiltonian
        self.tot_time=0          # orginal Hamiltonian
        self.Htarget=None
        self.randseed=0

        self.state=None
        self.frame=np.eye(2) # Toggling frame unitary
   

    def setParam(self,maxTime, nSpin, min_delay, max_delay, pw):
        self.maxTime=maxTime
        self.action_space=spaces.Discrete(45)# x y -x -y &  angle
        self.observation_space = spaces.Box(low=-1, high=1,shape=(13,), dtype=np.float32)
        self.nSpin=nSpin
        self.unitary=np.eye(2**nSpin) # Unitary operator
        self.pw=pw
        self.min_delay=min_delay
        self.max_delay=max_delay


    def setTargetH(self,target):
        if self.nSpin!=None:
            sh=target.shape
            if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
                raise ValueError('target must be a 2^nSpin-by-2^nSpin matrix.')
        else:
            raise ValueError('setParam before setTarget')
        self.targetH=target




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
        
    def setHtarget(self,Htarget):
        if self.nSpin!=None:
          sh=Htarget.shape
          if sh[0]!=2**self.nSpin or sh[1]!=2**self.nSpin:
            raise ValueError('U0 must be a 2^nSpin-by-2^nSpin matrix.')
        else:
          raise ValueError('setParam before setU0')
        self.Htarget=Htarget
        
        
        
        
    def Hbar(self,n):
        tauu=1
        e=1/np.pi/1.7
        h1=np.zeros([2**self.nSpin,2**self.nSpin])+1j*np.zeros([2**self.nSpin,2**self.nSpin])
        h2=0
        
        act=self.state[n]
        
        

        index=int(act/2)

        if act%2==0:
            phi=np.pi/2
        else:
            phi=np.pi/2+e

        field=self.H_pulse_nSpin(index)
            
        
        for i in range(2**self.nSpin):
            for j in range(2**self.nSpin):
                f = lambda x: np.array( np.dot(np.dot(expm(+1j*phi/2*x/self.pw*field),self.H0),expm(-1j*phi/2*x/self.pw*field))[i][j].real )
                g = lambda x: np.array( np.dot(np.dot(expm(+1j*phi/2*x/self.pw*field),self.H0),expm(-1j*phi/2*x/self.pw*field))[i][j].imag )
                h1[i][j]=integrate.quad(f, 0, self.pw)[0]+1j*integrate.quad(g, 0, self.pw)[0]
        
        u0=expm(-1j*phi/2*field)
        h2=np.dot(np.dot(np.transpose(np.conjugate(u0)),self.H0),u0)
        
        u=np.eye(2**self.nSpin)
        
        for pp in range(n-1):
            act=self.state[pp]
     
            index=int(act/2)
            
            if act%2==0:
                phi=np.pi/2
            else:
                phi=np.pi/2+e

            field=self.H_pulse_nSpin(index)
            
           
            u=np.dot(expm(-1j*phi/2*field),u)
        
        

        return np.dot(np.dot(np.transpose(np.conjugate(u)),h1+h2*(tauu)),u)/(tauu+self.pw)

    
    
    
    def getAHT1(self):
            H=0
            for p in range(6):
                H=H+self.Hbar(p)
            return H/6

    def getAHT2(self):
        H=0
        for p in range(6):
            for pp in range(p):
                H=H+commute(self.Hbar(p),self.Hbar(pp))
        return H*(self.pw+1)/2/1j/6

    def getAHT3(self):
        H=0
        for p in range(6):
            for pp in range(p):
                for ppp in range(pp):
                    H=H+commute(self.Hbar(p),commute(self.Hbar(pp),self.Hbar(ppp)))+commute(self.Hbar(ppp),commute(self.Hbar(pp),self.Hbar(p)))
        return -H*(self.pw+1)**2/6/6


    def getAHT4(self):
            H=0
            for p1 in range(6):
                for p2 in range(p1):
                    for p3 in range(p2):
                        for p4 in range(p3):
                            H=H+commute(commute(commute(self.Hbar(p1),self.Hbar(p2)),self.Hbar(p3)),self.Hbar(p4))+commute(self.Hbar(p1),commute(commute(self.Hbar(p2),self.Hbar(p3)),self.Hbar(p4)))+commute(self.Hbar(p1),commute(self.Hbar(p2),commute(self.Hbar(p3),self.Hbar(p4))))+commute(self.Hbar(p2),commute(self.Hbar(p3),commute(self.Hbar(p4),self.Hbar(p1))))
            return -H*(self.pw+1)**3/6/12/1j



    
    

# get Tr(self.unitary*self.target')
    def getFidelity(self):
        D, V=schur(self.unitary)
        dd=np.diag(np.exp(np.log(np.diag(D))/self.tot_time)) 
    #     dd=np.diag(np.exp(np.log(np.diag(D)))) 

        U=np.dot(np.dot(V,dd),np.transpose(np.conjugate(V))) # unitary operator in unit time
        return np.abs(np.trace(U*expm(1j*self.targetH))/2**self.nSpin)


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
    def H_pulse_nSpin(self, index):
        n=index
        if n==0:
            temp=Hamiltonian(self.nSpin,'p',[1],[1,0,0])
        elif n==1:
            temp=Hamiltonian(self.nSpin,'p',[1],[0,1,0])
        elif n==2:
            temp=Hamiltonian(self.nSpin,'p',[1],[-1,0,0])
        elif n==3:
            temp=Hamiltonian(self.nSpin,'p',[1],[0,-1,0])
        return temp
    
    
    
    def step(self, action,stepp,xy):
        pw=self.pw
        e=1/2/np.pi
        ee=0.2
        if action==44:
            field=np.eye(2**self.nSpin)
            phi=0
        else:

            index=int(action/11)
            dis=int(action%11)*ee
            phi=np.pi/2+dis*e
            field=self.H_pulse_nSpin(index)

            
        

#         phi=np.pi/2
        tauu=1
       
        self.tot_time=self.tot_time+tauu+pw*1
        
        self.state[stepp]=action/44
        self.state[12]=stepp/11
        
#         ##### angle imperfection
#         phi=phi+np.pi/300*(np.random.rand(1)[0]*2-1)
        phi=phi+np.pi/200
#         phi=phi+1/2/np.pi*0.4
#         phi=phi+np.pi/200*self.randseed
        
        
       
        if xy==0:
            unitary_n=np.dot(expm(-1j*self.H0*tauu),expm(-1j*(self.H0*pw+phi/2*field)))### pulse before eolution
        elif xy==1:
            unitary_n=np.dot(expm(-1j*(self.H0*pw+phi/2*field)),expm(-1j*self.H0*tauu))### pulse after evolution
        else:
            raise ValueError('xy must be 0 or 1')




        self.unitary=np.dot(unitary_n,self.unitary)
        
        if stepp>10:
            reward=-np.log(1-self.getFidelity())
            if self.tot_time<0.1:
                reward=-10
            return self.state, reward, True, {"frame":self.frame}

#         unitary_n=np.dot(np.dot(u_frame,self.U0),np.transpose(np.conjugate(u_frame)))
#         self.unitary=np.dot(unitary_n,self.unitary)

        return self.state, 0.00, False, {"frame":self.frame}      


  
    def reset(self):
        self.tot_time=0
        self.randseed=np.random.rand(1)[0]
        self.state=np.array([-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.])
        if self.nSpin==None:
            self.unitary=None # Unitary operator
        else:
            self.unitary=np.eye(2**self.nSpin)
        self.frame=np.eye(2)
        return self.state, {"frame":self.frame}

    def render(self, mode='human', close=False):
        pass
