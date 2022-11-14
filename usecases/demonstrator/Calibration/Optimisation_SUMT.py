# -----------------
# https://web.stanford.edu/class/ee364a/lectures/stoch_prog.pdf (good material for stochastic programming)
# Literature used for implementation:
# [1] : 1. Schulman, J., Heess, N., Weber, T. & Abbeel, P. Gradient Estimation Using Stochastic Computation Graphs. Preprint at http://arxiv.org/abs/1506.05254 (2016).
# [2] : 1. Wang, I.-J. & Spall, J. C. Stochastic optimization with inequality constraints using simultaneous perturbations and penalty functions. in 42nd IEEE International Conference on Decision and Control (IEEE Cat. No.03CH37475) 3808–3813 (IEEE, 2003). doi:10.1109/CDC.2003.1271742.
# [3] : 1. Bird, T., Kunze, J. & Barber, D. Stochastic Variational Optimization. Preprint at http://arxiv.org/abs/1809.04855 (2018).
# [4] : 1. Dimitriev, A. & Zhou, M. ARMS: Antithetic-REINFORCE-Multi-Sample Gradient for Binary Variables. Preprint at http://arxiv.org/abs/2105.14141 (2021).
# [5] : 1. Byrne, C. Sequential unconstrained minimization algorithms for constrained optimization. Inverse Problems 24, 015013 (2008).

import sys
sys.path.extend(['/home/atul/PhD_Tasks/LeBeDigital/ModelCalibration'])

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

# local imports
from usecases.demonstrator.StructuralSolver.Column_simulation import Column_simulation


# The script tries to just solve the optmization problem for column simulation.
# Doesnt target to modularize it yet.


# -- Load the calibrated parameters
# Note: Forgot to save the values, so just getting from the graph
#phi_mean = np.hstack((np.array([-0.7,0.045,0.009,-0.4]).reshape(-1,1), np.array([2.35, 6.25, 3.55, 4.24]).reshape(-1,1)))
#phi_sd = np.array([-7.5,-12.8,-11.4,-8.9])
phi_mean = np.load('usecases/demonstrator/Calibration/Results/phi_mean.npy')
phi_sd = np.load('usecases/demonstrator/Calibration/Results/phi_sd.npy')
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
    """
    Forward model interface. Inputs latens and outputs the KPIs need downstream for the optimization problem
    TODO: make it more general. Make it a class method which needs to be overloaded later.
    Args:
        b:

    Returns:

    """

    # (time, max temp for x = 0) =tensor([5673.6180, 80]), (time, max temp for x = 1.) = tensor([9198.5823,   56]). So choose temp value in between as constraint
    scaling = np.array([1e-04, 1e-03, 1, 1e05])
    latents = b.detach().numpy()*scaling
    data, time, temp = Column_simulation(latents)
    return th.as_tensor(np.array([time[0], temp])) # time is the point where yeild changes sign and temp is the max temp of the list

pr = Prior_(phi=phi_test)
#chk_1 = forward_model(pr._b_mean(th.tensor([0.]),th.from_numpy(phi_mean))) # tensor([5484.9441,   80.1733])
#chk_2 = forward_model(pr._b_mean(th.tensor([1.]),th.from_numpy(phi_mean))) # tensor([8905.9446,   56.8212])
#chk_3 = forward_model(pr._b_mean(th.tensor([0.5]),th.from_numpy(phi_mean))) # tensor([6909.0877,   68.6847])
#chk_4 = forward_model(pr._b_mean(th.tensor([0.6]),th.from_numpy(phi_mean)))
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

def objective(X: th.Tensor, C: th.Tensor):
    """
    Constructs the final "augemented" objective (converts constrained problem to an unconstrained one) to be passed /
    to an optimiser. Needs forward model, objective and constraints
    https://pytorch.org/docs/stable/distributions.html Talks about building a stochastic graph
    https://arxiv.org/pdf/1506.05254.pdf
    TODO: Make this a method to be overloaded later.
    TODO: Update for it to work with objective and constraints specific by function method not hardcoded.
    Args:
        X: The design/optimization variables
        C: The penalty term scaling factor

    Returns:
        object: tuple of augmented objective, objective, constraints and the stochastic values of the forward model
        output at x^*
    """
    assert isinstance(X,th.Tensor)
    assert X.requires_grad == True
    # Values which needs to be adjusted
    alpha = th.tensor(68) # The temp value which should not be exceeded. for x=0.5
    pr = Prior_(phi=phi_test)
    O_x = []
    C_x = []
    loss_aug =[]
    Y_b_N =[]
    N= 30 # no of samples for Monte Carlo estimates

    # -- Score function estimator
    # Monte carlo estimates
    for i in range(N): # E_{p(b|x,phi)} [y_o(b)]

        # -- sampling and calling forward solver for that sample
        b_sample = pr.sample(X,samples=1) # https://github.com/pytorch/pytorch/issues/7637 (>1 read this)
        val = pr.logeval(b_sample, x=X)
        forward_b = forward_model(b_sample.flatten())

        # Reinforce algo
        # E[(O(x) + c C(x)) log p(b|x,phi^star)]_p(b|x,phi^star)
        #TODO: below is hard coded and ugly. Make separate functions to do so.
        obj = forward_b[0]
        constraint_1 = forward_b[1]
        aug_obj_tmp = val*(obj + C*th.max(constraint_1-alpha,th.tensor(0)))
        #out = forward_b*val

        # data collection
        #prob_sum.append(val)
        Y_b_N.append(forward_b)
        O_x.append(obj) # passing time here
        C_x.append(constraint_1)
        loss_aug.append(aug_obj_tmp)
    #Z = th.sum(th.stack(prob_sum))
    #O_x_hat = th.sum(th.stack(O_x),axis=0)/Z
    #C_x_hat = th.sum(th.stack(C_x),axis=0)/Z

    # --- MC estimates
    aug_obj_hat = th.sum(th.stack(loss_aug))/N
    O_x_hat = th.sum(th.stack(O_x))/N
    C_x_hat = th.sum(th.stack(C_x))/N


    # -- Pathwise derivative (Works only when forward model is differentiable else no)
    # for i in range(N):
    #     forward_b = forward_model(b_samples[i, :])
    #     V_x.append(forward_b[0])
    #     C_x.append(forward_b[1])
    # V_x_hat = th.mean(th.stack(V_x))
    # C_x_hat = th.mean(th.stack(C_x))
    #obj =  O_x_hat + C*th.max(C_x_hat-alpha,th.tensor(0))
    #obj =coeff*val + alpha +b_sample
    assert aug_obj_hat.requires_grad == True
    return aug_obj_hat, O_x_hat, C_x_hat, Y_b_N

#X = th.tensor([0.8], requires_grad=True)
#tmp, a, b, c = objective(X)
#tmp.backward()
# print(X.grad)


def run(x_init:float,eps =0.001, eps_opt = 0.001, verbose = True) -> None:
    """
    https://mat.uab.cat/~alseda/MasterOpt/const_opt.pdf
    TODO: Make a method of the stochastic optimisazation class which inputs the "augmented" objective.
    Args:
        x_init: The initial value of the optimisation
        eps: The stopping crition for the SUMT algorigthm
        eps_opt: The stopping criterion for the inner loop optimization
        verbose:

    Returns:
        dataframe: A dataframe containing evolution of loss, objective, contraints, design variables,
        penalty scaling term
        Y_b: [ugly, need to update] The fowrad model output at x^*

    """
    X = th.tensor(x_init, requires_grad=True)
    C = th.tensor(5000,requires_grad=False)
    optimizer = th.optim.Adam([X], lr=0.01)
    losses = []
    objective_value = []
    constraints = []
    x_inmdt = [] # Intermediate for tracking
    grad = []
    c = [] # penalty parameter
    #Y_b_step = []
    num_steps = 100
    k = 0
    for j in range(10):
        X_tmp = X.clone().detach() # temp variable to be later ued in SUMT
        for i in range(num_steps):
            optimizer.zero_grad()
            # Y_b is the samples of the solver output for the last opt step.
            loss, O_x, C_x, Y_b = objective(X,C) # append with - sign if doing argmax
            #loss, O_x, C_x, Y_b = objective(X)
            loss.backward()
            optimizer.step()
            losses.append(loss)
            x_inmdt.append(X.clone())
            grad.append(X.grad.clone())
            objective_value.append(O_x)
            constraints.append(C_x)
            c.append(C)
            #Y_b_step.append(Y_b)

            if verbose:
                #if num_steps % 5 == 0:
                print(f"Iteration :{i+1}, loss value: {loss}, Objective : {O_x}, Constraints : {C_x}, x value: {X}, grad w.r.t x: {X.grad}, C : {C}")
            if i>0:
                if th.abs(X - x_inmdt[-2]) < eps_opt: # for covergance check of optimizer
                    print("----------------- The optimiser is Converged !! ----------------------")
                    break
        if th.abs(X_tmp - X)<eps:
            print(f"The SUMT is done. The final penalty paramter is {C}")
            break
        else:
            print(f"The SUMT is not done yet. The penalty paramter is {C} for the {k}th step")
            C = 1.5*C
            k+=1

            # else:
            #     C = 2*C


    data = {'loss':th.stack(losses).detach().numpy(),
            'X':th.cat(x_inmdt).detach().numpy(),
            'X_grad':th.cat(grad).detach().numpy(),
            'E_objective':th.stack(objective_value).detach().numpy(),
            'E_constraints': th.stack(constraints).detach().numpy(),
            'C': th.stack(c).detach().numpy()
            }
    df = pd.DataFrame(data=data)
    return df, Y_b

# sandboxing
if __name__ == "__main__":

    df, Y_b= run([0.2])

    df.to_csv('./OptimisationResults_'+datetime+'.csv')
    np.save('./Y_b_opt_x'+datetime+'.npy',th.stack(Y_b).detach().numpy())




