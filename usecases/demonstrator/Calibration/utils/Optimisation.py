# -----------------
# https://web.stanford.edu/class/ee364a/lectures/stoch_prog.pdf (good material for stochastic programming)
#


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import torch
from tqdm import tqdm

import torch as th
torch.set_default_dtype(torch.float64)

import os
from datetime import datetime

import matplotlib as mpl
from matplotlib import rc
mpl.rcParams['font.family'] = ['times new roman'] # default is sans-serif
rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

# -- Load the calibrated parameters



# -- Constrtust the prior on the latents p(b|x;\varphi)
class Prior_(object):
    def __init__(self, phi: list):

        self.phi = phi
        self.cov = None

    def _b_mean(self, x, WB):
        assert WB.ndim == 2
        b_vec = th.matmul(WB[:, :-1], x) + WB[:, -1]
        return b_vec

    def logeval(self, b, x):
        assert isinstance(x, th.Tensor)
        assert x.requires_grad == True
        phi_mean = self.phi[0]
        phi_sd_diag = self.phi[1]
        mean = self._b_mean(x, th.from_numpy(phi_mean))
        assert mean.shape[0] == phi_sd_diag.shape[0]
        phi_sd_diag = th.from_numpy(phi_sd_diag)  # diagonal entries of cov
        # self.cov = th.diag(phi_sd_diag_) @ th.diag(phi_sd_diag_).mT
        cov = th.diag(1e-07 + th.exp(phi_sd_diag))
        dist = th.distributions.MultivariateNormal(mean, cov)
        val = dist.log_prob(b)
        #val.backward()
        #grad_phi = phi_.grad
        #grad_sigma = phi_sd_diag_.grad
        # returing falttened gradients
        return val #, grad_phi,grad_sigma  # negative as later grad ascent needs to performed to find arg max logp(D|phi)

    def sample(self, x, samples=100):
        assert isinstance(x, th.Tensor)
        assert x.requires_grad == True

        phi_mean = self.phi[0]
        phi_sd_diag = self.phi[1]
        phi_sd_diag = th.from_numpy(phi_sd_diag)
        mean = self._b_mean(x, th.from_numpy(phi_mean))
        # cov = th.diag(phi_sd_diag) @ th.diag(phi_sd_diag).mT
        cov = th.diag(1e-07 + th.exp(phi_sd_diag))
        dist = th.distributions.MultivariateNormal(mean, cov)
        samples = dist.sample([samples, ])

        return samples


# -- Defind the structural model
def forward_model(b):

    temp = b.detach().numpy()
    # test function
    time = np.max(temp)
    temp = np.min(temp)
    return th.as_tensor([time, temp]) # time is the point where yeild changes sign and temp is the max temp of the list


# -- DEfining the optimisation problem

def V_x():
    """Define the obejctive here. Returns approximation of the expectaion."""
    return NotImplementedError

def C_x():
    """Define the contraints.Returns approximation of the expectaion."""
    return NotImplementedError

def MC_approx():
    """defining Monte Carlo approximation for the integrals. Use to to approximate the Expected objective
     and constraints"""
    return NotImplementedError

def objective(X):
    """Constructs the final objective to be passed to an optimiser with the V(x) and C(x)"""
    assert isinstance(X,th.Tensor)
    assert X.requires_grad == True
    # Values which needs to be adjusted
    alpha = th.tensor(65) # The temp value which should be exceeded
    coeff =1
    phi_mean = np.hstack((np.ones((4, 1)), np.array([2.916, 2.4229, 5.554, 5.0]).reshape(-1,1)))
    phi_sd = -1 * np.ones(4)
    phi_test = [phi_mean, phi_sd]
    pr = Prior_(phi=phi_test)
    V_x = []
    C_x = []
    prob_sum = []
    N= 100 # no of samples for Monte Carlo estimates
    b_samples = pr.sample(X,samples=N)
    # Monte carlo estimates
    for i in range(N): # E_{p(b|x,phi)} [y_o(b)]
        #b_sample = pr.sample(X,samples=1)
        #assert  b_sample.requires_grad == True
        val = th.exp(pr.logeval(b_samples[i,:],x=X)) # exp as it is logprob
        prob_sum.append(val)
        out = forward_model(b_samples[i,:])*val
        #print(X.grad)
        V_x.append(out[0]) # passing time here
        C_x.append(out[1])
    V_x_hat = th.sum(th.stack(V_x),axis=0)/th.sum(th.stack(prob_sum))
    C_x_hat = th.sum(th.stack(C_x),axis=0)/th.sum(th.stack(prob_sum))

    obj =  V_x_hat + coeff*th.min(C_x_hat,alpha)
    #obj =coeff*val + alpha +b_sample
    assert obj.requires_grad == True
    return obj

X = th.tensor([0.9], requires_grad=True)
tmp = objective(X)
# print(X.grad)


def run(x_init:float, verbose = True) -> None:
    X = th.tensor(x_init, requires_grad=True)
    optimizer = th.optim.Adam([X], lr=0.01)
    losses = []
    x_inmdt = [] # Intermediate for tracking
    grad = []
    num_steps = 40
    for i in range(num_steps):
        optimizer.zero_grad()
        loss = objective(X) # append with - sign if doing argmax
        loss.backward()
        # print(XX.grad)
        optimizer.step()
        losses.append(loss.item())
        x_inmdt.append(X)
        grad.append(th.norm(X.grad))

        if verbose:
            if num_steps % 10 == 0:
                print(f"Iteration :{i+1}, Objective value: {loss}, x value: {X}, grad w.r.t x: {X.grad} ")


# sandboxing
run([0.5])

# th.min(th.tensor(0.5),0.1)