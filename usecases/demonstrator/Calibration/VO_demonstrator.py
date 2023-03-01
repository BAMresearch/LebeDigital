# atul.agrawal@tum.de (Data Driven Materials Modeling Group)
# Trying to implement 1. Bird, T., Kunze, J. & Barber, D. Stochastic Variational Optimization.
# Preprint at http://arxiv.org/abs/1809.04855 (2018).(The Fig 2 specifically )
# 2. Other implementations for inspiration:
#  - https://github.com/aajanki/variational-optimization/blob/master/variationaloptimization/optimize.py
#  - https://github.com/artix41/AVO-pytorch/blob/master/avo-poisson.ipynb

# observartions/updates :
# 9.11.2022 : This code is working as it should. Accurately recreating the Fig 2 of the paper.
# 211.02.2023 : adding a constraints C(x) \geq 1 also. Check notes for derivation

import sys
import os

sys.path.extend(['/home/atul/PhD_Tasks/LeBeDigital/ModelCalibration'])  # temp fix to add the project path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch as th

th.set_default_dtype(th.float64)
import json
import os, sys
from datetime import datetime

import matplotlib as mpl
from matplotlib import rc

plt.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,bm}']
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

# local imports
from usecases.demonstrator.Calibration.utils.viz import plot_constraints_and_objective


def load_json(path: str) -> dict:
    if path[-5:] == '.json':
        with open(path) as f:
            data = json.load(f)
    return data


def update_json(file_path: str, key: str, value):
    # Read the JSON file
    with open(file_path, 'r') as f:
        data = json.load(f)
    # TODO:will work only when 'value' key is present
    # Update the value of the specified key
    data[key]['value'] = value

    # Write the updated data back to the JSON file
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4, sort_keys=True)


Optimization_workflow_path = '../../optimization_paper/optimization_workflow'
Results_path = Optimization_workflow_path + '/Results/'
# FEM_KPI = Results_path + 'kpi_from_fem.json'
# gwp_KPI = Results_path + 'gwp_beam.json'
# beam_design_KPI = Results_path + 'beam_design.json'

Input_path = Optimization_workflow_path + '/Inputs/'
aggregate_ratio_path = Input_path + 'aggregates_volume_fraction.json'
slag_ratio_path = Input_path + 'sc_fraction.json'  # sc slag/cement ratio. Instead of slag, it can be some type of cem too
phi_hydration_path = Input_path + 'phi_hydration.json'
phi_paste_path = Input_path + 'phi_paste.json'

X = {'agg_ratio': 0.6, 'slag_ratio': 0.4}
seed = 5


def function(X: dict, seed: int) -> dict:
    """
    Runs the snakemake workflow and the returns the KPIs for objective and constraints for a given value of the design
    variables. The Random variables (b) x->b->KPIs are also sampled for a given value of seed.
    Args:
     X: (dict) with keys 'agg_ratio' (volume ratio of the aggregates) and 'slag_ratio'
     seed: the seed parameter. This ensures that the sampled Random variable here is the same as the one passed in the
     forward call

    Returns:
        y : dict with all the KPIs

    """
    # Pass the parameter to X to the input to forward. Meaning overwrrite the input.
    # The design variables, aggregate ratio and the slag ratio needs to be updated.
    update_json(aggregate_ratio_path, 'aggregates_volume_fraction', X['agg_ratio'])
    update_json(slag_ratio_path, 'sc_volume_fraction', X['slag_ratio'])

    # pass the seed to the scripts for the RVs (see eqn 29 SVO paper)
    # Updating the phi's which are input to the script.
    update_json(phi_hydration_path, 'seed', seed)
    update_json(phi_paste_path, 'seed', seed)

    # Run the workflow using snakemake
    # add the path to the workflow file and the path to the directory
    workflow_file_path = Optimization_workflow_path + '/Snakefile'
    directory_path = Optimization_workflow_path
    os.system(f'snakemake --cores 6 --snakefile {workflow_file_path} '
              f'--directory {directory_path}  workflow_targets --use-conda')

    # Read in the KPIs in a dict
    FEM_KPI = Results_path + 'kpi_from_fem.json'
    gwp_KPI = Results_path + 'gwp_beam.json'
    beam_design_KPI = Results_path + 'beam_design.json'
    y = {}
    for i, path in enumerate([FEM_KPI, gwp_KPI, beam_design_KPI]):
        tmp = load_json(path)
        y.update(tmp)

    # return the KPIs
    return y

#tmp = function(X,seed)

class objective_constraints_demonstrator:
    def __init__(self, function: callable):
        self.function = function
        self.KPI_store = None

    def objective(self, x_1, x_2, **kwargs) -> float:
        seed = kwargs['seed']
        X_design= {'agg_ratio': x_1, 'slag_ratio': x_2}
        self.KPI_store = self.function(X_design, seed)
        y_o = self.KPI_store['gwp_mix']['value']  # the gwp of th beam
        return y_o

    def constraint(self, x_1, x_2, **kwargs) -> float:
        KPI = kwargs['KPI_yc']
        y_c = self.KPI_store[KPI]['value']
        return y_c

    def constraint_1(self, x_1, x_2, **kwargs) -> float:
        return self.constraint(x_1, x_2, KPI_yc='check_steel_area')

    def constraint_2(self, x_1, x_2, **kwargs) -> float:
        return self.constraint(x_1, x_2, KPI_yc='max_reached_temperature')

    def constraint_3(self, x_1, x_2, **kwargs) -> float:
        return self.constraint(x_1, x_2, KPI_yc='time_of_demolding')


# oc = objective_constraints(function)
# y_o = oc.objective(X,seed=3)
# y_c_1 = oc.constraint(KPI='check_steel_area')
# y_c_2 = oc.constraint(KPI='max_reached_temperature')
# y_c_3 = oc.constraint(KPI='time_of_demolding')
plot = False
if plot:
    x_bounds = (0.1, 0.9)
    y_bounds = (0.1, 0.9)
    oc = objective_constraints_demonstrator(function)
    constraints = [oc.constraint_1, oc.constraint_2, oc.constraint_3]
    fig_main, obj_val, cons_val = plot_constraints_and_objective(oc.objective, constraints, x_bounds, y_bounds,
                                                                 x_steps=5, y_steps=5, seed=2,
                                                                 KPI_yc='check_steel_area')

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))

    x_grid = np.linspace(x_bounds[0], x_bounds[1], 5)
    y_grid = np.linspace(y_bounds[0], y_bounds[1], 5)
    for k, con_vals_k in enumerate(cons_val):
        ax[k].set_aspect('equal')
        cons_contour = ax[k].contourf(
            x_grid,
            y_grid,
            con_vals_k,
            levels=10,
            # this ensures that just a sharp deviding line is present, kind of like indicator function. marks where there is a change of sign
            cmap='inferno'
        )
        ax[k].set_xlabel('$x_1$')
        ax[k].set_ylabel('$x_2$')
        ax[k].set_title('Constraint-' + str(k + 1))
        plt.colorbar(cons_contour)
    plt.show()

    fig.savefig('usecases/demonstrator/Calibration/Results/constraints_vs_design_variables' + datetime + '.pdf')
    fig.savefig('usecases/demonstrator/Calibration/Results/constraints_vs_design_variables' + datetime + '.png')
    fig, ax = plt.subplots()

    ax.set_aspect('equal')
    obj_contour = ax.contourf(x_grid, y_grid, obj_val, levels=20, cmap='inferno')
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.set_title('Objective Function')
    plt.colorbar(obj_contour)
    plt.show()
    fig.savefig('usecases/demonstrator/Calibration/Results/Objective_vs_design_variables' + datetime + '.pdf')
    fig.savefig('usecases/demonstrator/Calibration/Results/Objective_vs_design_variables' + datetime + '.png')


def MVN(mu: list, cov: list):
    # define the parametric mean
    dist = th.distributions.MultivariateNormal(th.as_tensor(mu), th.as_tensor(cov))
    return dist


# load \phi into dict

phi_hydration = load_json(phi_hydration_path)
phi_paste = load_json(phi_paste_path)


def objective(x_1, x_2, **kwargs):
    """
    # TODO: add a separate variational dist function
    Args:
        mu:
        sigma:
        beta:

    Returns:

    """
    if isinstance(x_2, dict):
        mean = x_2['mean']
        std = x_2['std']
        # dist for design variable wrt grad is not there
        q_x_2 = th.distributions.Normal(th.as_tensor(mean), th.as_tensor(std))

    # define dist of b_1
    # TODO: write a class/fn for relation between design variable and latents
    # get input from phi_hydration.json
    mu_b_1 = phi_hydration['phi_mean']['value']
    cov_b_1 = phi_hydration['phi_cov']['value']
    mean_b_1 = th.matmul(th.tensor(mu_b_1)[:, :-1], x_1) + th.tensor(mu_b_1)[:, -1]
    q_b_1 = MVN(mean_b_1, th.as_tensor(cov_b_1))

    # define dist of b_2
    # get input from phi_paste.json
    mu_b_2 = phi_paste['phi_mean']['value']
    cov_b_2 = phi_paste['phi_cov']['value']
    # mean_b_2 = th.matmul(th.tensor(mu_b_2), x_1)
    mean_b_2 = th.tensor(mu_b_2) * (x_1+th.tensor(0.5))
    q_b_2 = MVN(mean_b_2, th.as_tensor(cov_b_2))

    num_samples = kwargs['num_samples']

    # defining holders
    U_theta_holder = []

    for i in range(num_samples):
        # set seed
        random_seed = i+420
        th.manual_seed(random_seed)

        # collect RV samples
        b_1 = q_b_1.sample()
        b_2 = q_b_2.sample()
        if isinstance(x_2, dict):
            x_2 = q_x_2.sample()

        # intstance for objectives and constraints
        oc = objective_constraints_demonstrator(function)
        # define objetcive
        obj = oc.objective(x_1=x_1.item(), x_2=x_2.item(), seed=random_seed)

        # define constraints
        # --- Set inputs for the constraints
        time_max = th.tensor(3)
        temp_max = th.tensor(70)
        c_1 = 1
        c_2 = 1
        c_3 = 1
        C_x_1 = oc.constraint_1(x_1, x_2)
        G_x_1 = c_1 * th.max(-th.as_tensor(C_x_1), th.tensor(0))
        C_x_2 = oc.constraint_2(x_1, x_2)
        G_x_2 = c_2 * th.max(th.as_tensor(C_x_2) - temp_max, th.tensor(0))
        C_x_3 = oc.constraint_3(x_1, x_2)
        G_x_3 = c_3 * th.max(th.as_tensor(C_x_3) - time_max, th.tensor(0))
        constraints = G_x_1 + G_x_2 + G_x_3

        # with constraints
        U_theta_holder.append((obj + constraints) * (q_b_1.log_prob(b_1) + q_b_2.log_prob(b_2) + q_x_2.log_prob(x_2)))
        # w/o constraints
        # U_theta_holder.append((th.as_tensor(y)) * dist.log_prob(x_sample))
    U_theta = th.sum(th.stack(U_theta_holder)) / num_samples

    assert U_theta.requires_grad == True
    return U_theta

# check
tmp = objective(x_1=th.tensor([0.4], requires_grad=True), x_2={'mean': [0.4], 'std': [0.1]}, num_samples=2)

def optimize(mu_init: float, eps=0.001, verbose=True) -> None:
    mu = th.tensor(mu_init, requires_grad=True)
    sigma = th.tensor([5.])
    beta = th.tensor(2 * th.log(sigma), requires_grad=True)
    # C = th.tensor(50,requires_grad=False)
    optimizer = th.optim.SGD([mu, beta], lr=0.1)
    losses = []
    objective_value = []
    constraints = []
    x_inmdt = []  # Intermediate for tracking
    sigma_list = []
    grad = []
    # Y_b_step = []
    num_steps = 150
    for i in range(num_steps):
        optimizer.zero_grad()
        # Y_b is the samples of the solver output for the last opt step.
        # loss, O_x, C_x, Y_b = objective(X,C) # append with - sign if doing argmax
        loss = objective(mu, sigma, beta=beta)
        # compute grads
        loss.backward()
        # print(XX.grad)
        losses.append(loss)
        x_inmdt.append(mu.clone())
        # sigma_list.append(sigma)
        sigma_list.append(th.sqrt(th.exp(beta.clone())))
        grad.append(th.norm(mu.grad.clone()))
        optimizer.step()

        # Y_b_step.append(Y_b)

        if verbose:
            # if num_steps % 5 == 0:
            print(
                f"Iteration :{i + 1}, loss value: {loss}, mu value: {mu}, sigma value: {sigma},grad w.r.t x: {mu.grad} ")
        if i > 0:
            if th.norm(mu - x_inmdt[-2]) < eps:
                print("----------------- Converged !! ----------------------")
                break
    # data = {'loss':th.stack(losses).detach().numpy(),
    #         'X':th.cat(x_inmdt).detach().numpy(),
    #         'X_grad':th.stack(grad).detach().numpy(),
    #         }
    # df = pd.DataFrame(data=data)
    return th.stack(x_inmdt).detach().numpy(), th.stack(sigma_list).detach().numpy()


mu_evolution_1, sigma_evolution_1 = optimize(mu_init=[4., -4.])
mu_evolution_2, sigma_evolution_2 = optimize(
    mu_init=[-4., 0.])  # starting from constraint violation and crossing the optima

x = np.arange(-5.0, 5.0, 0.1)
y = np.arange(-5.0, 5.0, 0.1)
X, Y = np.meshgrid(x, y)  # grid of point
Z = function(X, Y)  # evaluation of the function on the grid

fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)


def plot_evolution(mu, sigma, color, fig, ax):
    ax[0].contourf(X, Y, Z, levels=20)
    ax[0].plot(mu[:, 0], mu[:, 1], 'x', color=color)
    ax[0].set_xlabel('$x_1$')
    ax[0].set_ylabel('$x_2$')
    ax[1].plot(sigma)
    ax[1].set_ylabel('$\sigma$')
    ax[1].set_xlabel('iterations')
    # plt.savefig('./Figs/theta_evolution_VO_' + datetime + '.pdf')
    plt.show()
    return fig


ax[0].contourf(X, Y, Z, levels=20)
ax[0].plot(mu_evolution_1[:, 0], mu_evolution_1[:, 1], 'x', color='r')
ax[0].plot(mu_evolution_2[:, 0], mu_evolution_2[:, 1], 'x', color='y')
ax[0].set_xlabel('$x_1$')
ax[0].set_ylabel('$x_2$')
ax[1].plot(sigma_evolution_1, 'r')
ax[1].plot(sigma_evolution_2, 'y')
ax[1].set_ylabel('$\sigma$')
ax[1].set_xlabel('iterations')
plt.savefig('./Figs/theta_evolution_VO_constraints_' + datetime + '.pdf')
plt.show()

plot_evolution(mu_evolution_1, sigma_evolution_1, 'r', fig, ax)
plot_evolution(mu_evolution_2, sigma_evolution_2, 'g', fig, ax)

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
