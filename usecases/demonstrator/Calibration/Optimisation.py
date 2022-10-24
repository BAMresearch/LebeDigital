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

# local imports
from usecases.demonstrator.StructuralSolver.Column_simulation import Column_simulation


# -- Load the calibrated parameters
# Note: Forgot to save the values, so just getting from the graph
phi_mean = np.hstack((np.array([-0.7,0.045,0.009,-0.4]).reshape(-1,1), np.array([2.35, 6.25, 3.55, 4.24]).reshape(-1,1)))
phi_sd = np.array([-3.5,-3.8,-3.4,-3.8])
phi_test = [phi_mean, phi_sd]


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
        #samples = dist.rsample([samples, ]) # If reparam can be done, needs differentiable solver
        return samples


# -- Defind the structural model
def forward_model(b):

    # temp = b.detach().numpy()
    # # test function
    # time = np.max(temp)
    # temp = np.min(temp)

    # (time, max temp for x = 0) =tensor([5673.6180, 80]), (time, max temp for x = 1.) = tensor([9198.5823,   56]). So choose temp value in between as constraint
    scaling = np.array([1e-04, 1e-03, 1, 1e05])
    latents = b.detach().numpy()*scaling
    data, time, temp = Column_simulation(latents)
    return th.as_tensor(np.array([time[0], temp])) # time is the point where yeild changes sign and temp is the max temp of the list

pr = Prior_(phi=phi_test)
#chk_1 = forward_model(pr._b_mean(th.tensor([0.]),th.from_numpy(phi_mean))) # tensor([5484.9441,   80.1733])
#chk_2 = forward_model(pr._b_mean(th.tensor([1.]),th.from_numpy(phi_mean))) # tensor([8905.9446,   56.8212])
#chk_3 = forward_model(pr._b_mean(th.tensor([0.5]),th.from_numpy(phi_mean))) # tensor([6909.0877,   68.6847])
chk_4 = forward_model(pr._b_mean(th.tensor([0.6]),th.from_numpy(phi_mean)))
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
    """Constructs the final objective to be passed to an optimiser with the V(x) and C(x)
    https://pytorch.org/docs/stable/distributions.html Talks about building a stochastic graph
    https://arxiv.org/pdf/1506.05254.pdf
    """
    assert isinstance(X,th.Tensor)
    assert X.requires_grad == True
    # Values which needs to be adjusted
    alpha = th.tensor(68) # The temp value which should not be exceeded. for x=0.5
    coeff =th.tensor(100)
    # phi_mean = np.hstack((np.ones((4, 1)), np.array([2.916, 2.4229, 5.554, 5.0]).reshape(-1,1)))
    # phi_sd = -1 * np.ones(4)
    # phi_test = [phi_mean, phi_sd]
    pr = Prior_(phi=phi_test)
    V_x = []
    C_x = []
    prob_sum = []
    tmp =[]
    N= 10 # no of samples for Monte Carlo estimates
    b_samples = pr.sample(X,samples=N)
    # -- Score function estimator
    # Monte carlo estimates
    for i in range(N): # E_{p(b|x,phi)} [y_o(b)]
        val = th.exp(pr.logeval(b_samples[i,:],x=X)) # exp as it is logprob
        #val = pr.logeval(b_samples[i, :], x=X)
        prob_sum.append(val)
        forward_b = forward_model(b_samples[i,:])
        out = forward_b*val
        tmp.append(forward_b)
        #print(X.grad)
        V_x.append(out[0]) # passing time here
        C_x.append(out[1])
    Z = th.sum(th.stack(prob_sum))
    V_x_hat = th.sum(th.stack(V_x),axis=0)/Z
    C_x_hat = th.sum(th.stack(C_x),axis=0)/Z

    # -- Pathwise derivative (Works only when forward model is differentiable else no)
    # for i in range(N):
    #     forward_b = forward_model(b_samples[i, :])
    #     V_x.append(forward_b[0])
    #     C_x.append(forward_b[1])
    # V_x_hat = th.mean(th.stack(V_x))
    # C_x_hat = th.mean(th.stack(C_x))
    obj =  0.1*V_x_hat + coeff*th.max(C_x_hat-alpha,th.tensor(0))
    #obj =coeff*val + alpha +b_sample
    assert obj.requires_grad == True
    return obj

# X = th.tensor([0.1], requires_grad=True)
# tmp = objective(X)
# print(X.grad)


def run(x_init:float,eps =0.005, verbose = True) -> None:
    X = th.tensor(x_init, requires_grad=True)
    optimizer = th.optim.Adam([X], lr=0.05)
    losses = []
    x_inmdt = [] # Intermediate for tracking
    grad = []
    num_steps = 18
    for i in range(num_steps):
        optimizer.zero_grad()
        loss = objective(X) # append with - sign if doing argmax
        loss.backward()
        # print(XX.grad)
        optimizer.step()
        losses.append(loss.item())
        x_inmdt.append(X.copy())
        grad.append(th.norm(X.grad))

        if verbose:
            #if num_steps % 5 == 0:
            print(f"Iteration :{i+1}, Objective value: {loss}, x value: {X}, grad w.r.t x: {X.grad} ")
        # if np.abs(X - x_inmdt[-2]) < eps:
        #     print("----------------- Converged !! ----------------------")
        #     break
    return losses, x_inmdt, grad

# sandboxing
loss, x, grad = run([0.8])

plt.plot(grad)
plt.plot(x)
plt.plot(loss)

# th.min(th.tensor(0.5),0.1)

import pandas

df = pandas.read_csv()