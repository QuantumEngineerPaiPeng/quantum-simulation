### Created by Xiaoyang Huang Jun,1st

import torch
import numpy as np

class  grad_inverter(object):
    """docstring for  grad_inverter:"""
    def __init__(self, action_bounds, grad, action):
        super( grad_inverter, self).__init__()
        self.action_size = len(action_bounds[0])
        action_bounds = np.asarray(action_bounds)
        self.action_input = action   
        self.pmax = torch.from_numpy(action_bounds[0]).type(torch.FloatTensor)
        self.pmin = torch.from_numpy(action_bounds[1]).type(torch.FloatTensor)
        self.prange = torch.from_numpy(np.asarray([x - y for x, y in zip(action_bounds[0],action_bounds[1])])).type(torch.FloatTensor)
        self.pdiff_max = torch.div(-self.action_input+self.pmax,self.prange)
        self.pdiff_min = torch.div(-self.action_input-self.pmin,self.prange)
        self.zeros_act_grad_filter = torch.from_numpy(np.zeros(self.action_size)).type(torch.FloatTensor)
        
        self.act_grad = grad[0]
        _, self.ones = torch.cat((self.act_grad, self.zeros_act_grad_filter), dim = -1).max(-1)
        self.grad_inverter_xy = self.ones * self.act_grad * self.pdiff_max + self.ones * self.act_grad * self.pdiff_min
