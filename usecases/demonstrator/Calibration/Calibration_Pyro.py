# Trying to get the phi and bs by pyro. Helpful links

# 1. http://pyro.ai/examples/svi_part_ii.html
# - The above has good intro how to wrap local (bs for diff datset N's) and global latents (phi's)
# - Nice discussion about conditional independance with pyro
# 2.  https://pyro.ai/examples/boosting_bbvi.html
# - If BBVI needs to be used, this serves as an example
# 3. https://forum.pyro.ai/t/sampling-multiple-observations-from-likelihood/4027
# - using MVN and plating example
# Note : If things dont work, put psuedocode and questions in the pyro forum

import yaml
import math
import os
import torch
import torch.distributions.constraints as constraints
import pyro
from pyro.optim import Adam
from pyro.infer import SVI, Trace_ELBO, TraceGraph_ELBO
import pyro.distributions as dist
from pyro.infer.autoguide import AutoDiagonalNormal

import fenics_concrete

# Generate observed data/ store exp data
data_file = '../artificial_hydration_data/artificial_hydration_data.yaml'
#Example 1:
# read file and access artificial data:


def data_for_inference(data_file):
    """takes the hydration data and genetares a form which can be used"""
    with open(data_file) as file:
        hydration_data = yaml.safe_load(file)
    return hydration_data

# Define forward model

def forward_model():
    parameter = fenics_concrete.Parameters()  # using the current default values

    # -- latents -----
    # parameter['B1'] = 2.916E-4  # in 1/s (le 0, < 0.1)
    # parameter['B2'] = 0.0024229  # - (le 0, smaller 1)
    # parameter['eta'] = 5.554  # something about diffusion (should be larger 0)
    # parameter['T_ref'] = 25  # reference temperature in degree celsius
    # parameter['Q_pot'] = 500e3 # potential heat per weight of binder in J/kg

    # -- adding scaling back the values
    parameter['B1'] = inp_latents[0] * 1e-04  # in 1/s (le 0, < 0.1)
    parameter['B2'] = inp_latents[1] * 1e-03  # - (le 0, smaller 1)
    parameter['eta'] = inp_latents[2]  # something about diffusion (should be larger 0)
    parameter['Q_pot'] = inp_latents[3] * 1e05  # potential heat per weight of binder in J/kg

    # -- observed inputs
    parameter['igc'] = 8.3145  # ideal gas constant in [J/K/mol], CONSTANT!!!
    parameter['zero_C'] = 273.15  # in Kelvin, CONSTANT!!!
    parameter[
        'E_act'] = 47002  # activation energy in Jmol^-1 (no relevant limits) (Depends only on simulated temp, if that is not change no need to infer E_act)
    parameter['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c (larger 0 and max 1)
    parameter['T_ref'] = 25  # reference temperature in degree celsius

    # this is the minimal time step used in the simulation
    # using a larger value will increase the speed but decrease the accuracy
    dt = 300  # value in seconds

    # this is the simulated temperature, needs to be adjusted depending on the temperature of the experimental data
    T = inp_obs['T_rxn']  # can be 20,40,60 as pert the exp values
    # this is the list of measured time data as given by the experiments
    # time_list = [0,5000,10000,20000,100000]
    time_list = inp_obs['time_list']

    # initiate material problem, for this the "fenics_concrete" conda package needs to be installed
    # use: 'mamba install -c etamsen fenics_concrete"
    problem = fenics_concrete.ConcreteThermoMechanical()

    # get the hydration function
    # this might change in the future to make it more easily accessible but for now it should work like this
    hydration_fkt = problem.get_heat_of_hydration_ftk()
    # the results are a heat list and a degree of hydration list, which you can ignore for now
    heat_list, doh_list = hydration_fkt(T, time_list, dt, parameter)

    return heat_list

chk = torch.diag(0.2*torch.tensor([2.916E-4, 0.0024229, 5.554, 500e3]))
# define probabilistic model
# TODO: Think how to pass the data, can be split into, x_hat, (time_step, y_hat: each with size Nx timesteps
def model(time_list, y_hat):
    # define global variable phi

    # --- if phi doesnt exist and we do p(b,data) = p(data|b) p(b)
    mean = torch.tensor([2.916, 0.0024229, 5.554, 500])
    cov = torch.diag(0.2 * mean)
    b = pyro.sample("\bm{b}", dist.MultivariateNormal(mean, covariance_matrix=cov))
    # define plate context (https://docs.pyro.ai/en/1.8.2/primitives.html), can be vectorised or serialized
    with pyro.plate("data",y_hat.shape[0]): # data can be N x timestep with N=5 here.
        # define MVN dist sample site for b,s --------------------------------------

        # --- if phi exists and we do p(b,phi,data) = p(data|b) p(b|phi)p(phi)
        # mean =
        # cov =
        # b = pyro.sample("\bm{b}",dist.MultivariateNormal(mean,covariance_matrix=cov))

        # call the solver with the bs -------------------------------------------------
        # Note if it throws differentiability error, use a offline trained surrogate here (maybe check BBVI too).
        y = forward_model()
        cov = torch.diag(1e-04*y)# define as much confidance on data, to start with can be 1% error
        # define the likelihood --------------------------------------------------------
        pyro.sample("\hat{y}",dist.MultivariateNormal(y,cov), obs=y_hat)



# vizualize the model to check
pyro.render_model(model, model_args=(data,), filename='./probabilistic_graph.pdf')

# define the variational dist for all the latents
def guide(data):

# -- or to simply things, use autoguide
guide = AutoDiagonalNormal(model)

# setup the optimizer
adam_params = {"lr": 0.0005, "betas": (0.90, 0.999)}
optimizer = Adam(adam_params)

# setup the inference algorithm
svi = SVI(model, guide, optimizer, loss=TraceGraph_ELBO())

# do gradient steps
pyro.clear_param_store()
for step in range(n_steps):
    loss = svi.step(data) # pass the data in appropraiet format, the same should be and arg for model and guide
    if step % 100 == 0:
        print("[iteration %04d] loss: %.4f" % (i + 1, loss / len(data)))
        print('.', end='')
        # TODO: add more diagnostics like gradients

# grab the learned variational parameters
para_1 = pyro.param("para_1").item()
para_2 = pyro.param("para_2").item()