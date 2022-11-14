# Trying to implement 1. Bird, T., Kunze, J. & Barber, D. Stochastic Variational Optimization.
# Preprint at http://arxiv.org/abs/1809.04855 (2018).(The Fig 2 specifically )
# 2. Other implementations for inspiration:
#  - https://github.com/aajanki/variational-optimization/blob/master/variationaloptimization/optimize.py
#  - https://github.com/artix41/AVO-pytorch/blob/master/avo-poisson.ipynb
# observartion : 9.11.2022 : This code is working as it should. Accurately recreating the Fig 2 of the paper.
import sys
sys.path.extend(['/home/atul/PhD_Tasks/LeBeDigital/ModelCalibration']) # temp fix to add the project path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from tqdm import tqdm

import torch as th
th.set_default_dtype(th.float64)

import os
from datetime import datetime

import matplotlib as mpl
from matplotlib import rc
mpl.rcParams['font.family'] = ['times new roman'] # default is sans-serif
rc('text', usetex=False)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

def function(X,Y):
    """

    Args:
     x:

    Returns:

    """
    # The function is not differentiable,
    # for D dimentional quadratic function
    #assert x.ndim ==1

    # Eq 49 in paper for 2d
    #y = 1/(2*len(x))*(np.sum([x[i]**2 for i in range(len(x))]))
    y = 1/(2*2)*(X**2 +Y**2)
    return np.array(y)

# y = function(2,1)
# plotting the function
x = np.arange(-5.0,5.0,0.1)
y = np.arange(-5.0,5.0,0.1)
X,Y = np.meshgrid(x, y) # grid of point
Z = function(X, Y) # evaluation of the function on the grid

im = plt.contourf(X,Y,Z, levels =20)
plt.show()

# class plots:
#     @staticmethod
    # line up plots here

def objective(mu,sigma, beta = None):
    """
    # TODO: add a separate variational dist function
    Args:
        mu:
        sigma:
        beta:

    Returns:

    """
    #assert theta.requires_grad == True
    # defining the va dist here
    #mu = theta[:-1] # \mu
    #sigma = theta[-1] # \sigma^2
    if beta is not None:
        #beta = 2*th.log(sigma)
        dist = th.distributions.MultivariateNormal(mu,th.exp(beta)*th.eye(2))
    else:
        dist = th.distributions.MultivariateNormal(mu, sigma**2 * th.eye(2))

    num_samples = 10
    #obj = th.zeros(num_samples)
    U_theta_holder = []
    for _ in range(num_samples):
        x_sample = dist.sample()
        y = function(x_sample[0],x_sample[1])
        U_theta_holder.append(th.as_tensor(y)*dist.log_prob(x_sample))
    U_theta = th.sum(th.stack(U_theta_holder))/num_samples

    assert U_theta.requires_grad == True
    return U_theta

#theta_check = th.tensor([4,-4, 5.], requires_grad=True)
#mu = th.tensor([0.,0.], requires_grad=True)
#sigma = th.tensor([1.])
#U_theta = objective(mu,sigma)

def optimize(mu_init:float,eps =0.001, verbose = True) -> None:
    mu = th.tensor(mu_init, requires_grad=True)
    sigma = th.tensor([5.])
    beta = th.tensor(2 * th.log(sigma),requires_grad=True)
    #C = th.tensor(50,requires_grad=False)
    optimizer = th.optim.SGD([mu,beta], lr=0.1)
    losses = []
    objective_value = []
    constraints = []
    x_inmdt = [] # Intermediate for tracking
    sigma_list = []
    grad = []
    #Y_b_step = []
    num_steps = 150
    for i in range(num_steps):
        optimizer.zero_grad()
        # Y_b is the samples of the solver output for the last opt step.
        #loss, O_x, C_x, Y_b = objective(X,C) # append with - sign if doing argmax
        loss = objective(mu,sigma,beta=beta)
        # compute grads
        loss.backward()
        # print(XX.grad)
        losses.append(loss)
        x_inmdt.append(mu.clone())
        #sigma_list.append(sigma)
        sigma_list.append(th.sqrt(th.exp(beta.clone())))
        grad.append(th.norm(mu.grad.clone()))
        optimizer.step()

        #Y_b_step.append(Y_b)

        if verbose:
            #if num_steps % 5 == 0:
            print(f"Iteration :{i+1}, loss value: {loss}, mu value: {mu}, sigma value: {sigma},grad w.r.t x: {mu.grad} ")
        if i>0:
            if th.norm(mu - x_inmdt[-2]) < eps:
                print("----------------- Converged !! ----------------------")
                break
    # data = {'loss':th.stack(losses).detach().numpy(),
    #         'X':th.cat(x_inmdt).detach().numpy(),
    #         'X_grad':th.stack(grad).detach().numpy(),
    #         }
    # df = pd.DataFrame(data=data)
    return th.stack(x_inmdt).detach().numpy(), th.stack(sigma_list).detach().numpy()



mu_evolution, sigma_evolution = optimize(mu_init=[4.,-4.])

x = np.arange(-5.0,5.0,0.1)
y = np.arange(-5.0,5.0,0.1)
X,Y = np.meshgrid(x, y) # grid of point
Z = function(X, Y) # evaluation of the function on the grid

im = plt.contourf(X,Y,Z, levels =20)
plt.plot(mu_evolution[:,0],mu_evolution[:,1],'x',color='r' )
plt.show()

plt.plot(sigma_evolution)
plt.show()

# class VO:
#     def __init__(self):
#
#     def objective:
#
#     def var_dist:
#
#
#     def run(self):


# --