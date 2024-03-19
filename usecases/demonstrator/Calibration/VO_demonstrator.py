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
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,bm}']
datetime = datetime.now().strftime("%d_%m_%Y-%I_%M_%S_%p")

# local imports
from usecases.demonstrator.Calibration.utils.viz import plot_constraints_and_objective
from lebedigital.demonstrator_optimization_scripts.utils import python_fn_run_jobs, read_kpis

# # %%
# def load_json(path: str) -> dict:
#     if path[-5:] == '.json':
#         with open(path) as f:
#             data = json.load(f)
#     return data


# # %%
# def update_json(file_path: str, key: str, value):
#     # Read the JSON file
#     with open(file_path, 'r') as f:
#         data = json.load(f)
#     # TODO:will work only when 'value' key is present
#     # Update the value of the specified key
#     data[key]['value'] = value

#     # Write the updated data back to the JSON file
#     with open(file_path, 'w') as f:
#         json.dump(data, f, indent=4, sort_keys=True)


# Optimization_workflow_path = '../../optimization_paper/optimization_workflow'
# Results_path = Optimization_workflow_path + '/Results/'
# # FEM_KPI = Results_path + 'kpi_from_fem.json'
# # gwp_KPI = Results_path + 'gwp_beam.json'
# # beam_design_KPI = Results_path + 'beam_design.json'

# Input_path = Optimization_workflow_path + '/Inputs/'
# aggregate_ratio_path = Input_path + 'aggregates_volume_fraction.json'
# slag_ratio_path = Input_path + 'sc_fraction.json'  # sc slag/cement ratio. Instead of slag, it can be some type of cem too
# phi_hydration_path = Input_path + 'phi_hydration.json'
# phi_paste_path = Input_path + 'phi_paste.json'

# X = {'agg_ratio': 0.6, 'slag_ratio': 0.4}
# seed = 5
# tmp = function(X,seed)

class objective_constraints_demonstrator:
    def __init__(self, function: callable):
        self.function = function
        self.KPI_store = None

    def objective(self, x_1, x_2, **kwargs) -> float:
        seed = kwargs['seed']
        X_design = {'slag_ratio': x_1, 'agg_ratio': x_2}
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
                                                                 x_steps=5, y_steps=5, seed=2)

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))

    x_grid = np.linspace(x_bounds[0], x_bounds[1], 5)
    y_grid = np.linspace(y_bounds[0], y_bounds[1], 5)
    X, Y = np.meshgrid(x_grid, y_grid)
    for k, con_vals_k in enumerate(cons_val):
        ax[k].set_aspect('equal')
        cons_contour = ax[k].contourf(
            X,
            Y,
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

    fig.savefig('./Results/constraints_vs_design_variables' + datetime + '.pdf')
    fig.savefig('./Results/constraints_vs_design_variables' + datetime + '.png')
    fig, ax = plt.subplots()

    ax.set_aspect('equal')
    obj_contour = ax.contourf(X, Y, obj_val, levels=20, cmap='inferno')
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.set_title('Objective Function')
    plt.colorbar(obj_contour)
    plt.show()
    fig.savefig('./Results/Objective_vs_design_variables' + datetime + '.pdf')
    fig.savefig('./Results/Objective_vs_design_variables' + datetime + '.png')

optimization = True
if optimization:
    def MVN(mu: list, cov: list):
        # define the parametric mean
        dist = th.distributions.MultivariateNormal(th.as_tensor(mu), th.as_tensor(cov))
        return dist


    # load \phi into dict

    #phi_hydration = load_json(phi_hydration_path)
    #phi_paste = load_json(phi_paste_path)

    def _translate_design_variable_to_stochastic(x:dict):
        """

        Parameters
        ----------
        x (dict) with keys:
        "mean" : the mean of the dist, assumed to be normal
        "s.d" : the std of the dist

        Returns
        -------
        q_x :  th.dist object
        """
        mean = x['mean']
        sd = x['s.d']
        assert mean.requires_grad == True, "The computational graph seems to be detached"
        assert sd.requires_grad == True, "The computational graph seems to be detached"
        q_x = th.distributions.Normal(loc=mean,scale=sd) #TODO: it just assumes normal now, can be log normal too.
        return q_x

    def _p_b_given_x(phi, x):
        """
        constract the stochstic node for b's for given phi and design variable x (slag)
        Parameters
        ----------
        phi: (dict) should contain the keys "phi_mean""value", "phi_cov""value" (its nested dict)
        x: the design variable (x_1 here)

        Returns
        -------

        """
        mu_b = phi['phi_mean']['value']
        cov_b = phi['phi_cov']['value']
        mean_b = th.matmul(th.tensor(mu_b)[:, :-1], x) + th.tensor(mu_b)[:, -1]
        q_b = MVN(mean_b,th.as_tensor(cov_b))
        return q_b


    # def objective(x_1, x_2, **kwargs):
        """
        # TODO: add a separate variational dist function
        Args:
            mu:
            sigma:
            beta:

        Returns:

        """
        if isinstance(x_2, dict):
            # mean = x_2['mean']
            # std = (x_2['std'])  # sigma is the input not sigma^2
            # # dist for design variable wrt grad is not there
            # assert mean.requires_grad == True
            # q_x_2 = th.distributions.Normal(mean, std)  # can use lognormal maybe to ensure +
            # print("using log Normal for x_2")
            q_x_2 = _translate_design_variable_to_stochastic(x=x_2)
            # q_x_2 = th.distributions.LogNormal(th.log(mean),std)
            # q_x_2 = th.distributions.LogNormal(mean, std)

        if isinstance(x_1,dict):
            q_x_1 = _translate_design_variable_to_stochastic(x=x_1)


        # define dist of b_1
        # get input from phi_hydration.json
        # mu_b_1 = phi_hydration['phi_mean']['value']
        # cov_b_1 = phi_hydration['phi_cov']['value']
        # mean_b_1 = th.matmul(th.tensor(mu_b_1)[:, :-1], x_1) + th.tensor(mu_b_1)[:, -1]
        # q_b_1 = MVN(mean_b_1, th.as_tensor(cov_b_1))
        #
        # # define dist of b_2
        # # get input from phi_paste.json
        # mu_b_2 = phi_paste['phi_mean']['value']
        # cov_b_2 = phi_paste['phi_cov']['value']
        # # mean_b_2 = th.matmul(th.tensor(mu_b_2), x_1)
        # # TODO: dummy now, need to write a proper function
        # # mean_b_2 = th.tensor(mu_b_2) * (x_1+th.tensor(0.5))
        # mean_b_2 = th.matmul(th.tensor(mu_b_2)[:, :-1], x_1) + th.tensor(mu_b_2)[:, -1]
        # q_b_2 = MVN(mean_b_2, th.as_tensor(cov_b_2))

        num_samples = kwargs['num_samples']

        # defining holders
        U_theta_holder = []
        obj_holder = []
        C_1_holder = []
        C_2_holder = []
        C_3_holder = []
        C_4_holder = []

        for i in range(num_samples):
            # set seed
            # The seed will ensure that the same RV samples are passed inside the forward model
            random_seed = np.random.randint(666)
            th.manual_seed(random_seed)

            # collect RV samples

            x_2 = q_x_2.sample()

            x_1 = q_x_1.sample()
            q_b_1 = _p_b_given_x(phi=phi_hydration,x=x_1)
            q_b_2 = _p_b_given_x(phi=phi_paste,x=x_1)
            b_1 = q_b_1.sample()
            b_2 = q_b_2.sample()

            # intstance for objectives and constraints
            oc = objective_constraints_demonstrator(function)
            # define objetcive
            # logistic sigmoid function to bound the input in 0-1 = 1/(1+e^(-y))
            x_1_scaled = th.special.expit(x_1)
            x_2_scaled = th.special.expit(x_2)
            obj = oc.objective(x_1=x_1_scaled.item(), x_2=x_2_scaled.item(), seed=random_seed)

            # define constraints
            # --- Set inputs for the constraints
            time_max = th.tensor(3)
            temp_max = th.tensor(70)
            max_agg_ratio = th.tensor(0.7)
            # workability constraint. Now temp that agg ratio < 0.6
            c_1 = 1e03
            c_2 = 0.1
            c_3 = 1
            c_4 = 1
            C_x_1 = oc.constraint_1(x_1_scaled, x_2_scaled)  # design criterion
            G_x_1 = c_1 * th.max(-th.as_tensor(C_x_1), th.tensor(0))
            C_x_2 = oc.constraint_2(x_1_scaled, x_2_scaled)  # temp
            G_x_2 = c_2 * th.max(th.as_tensor(C_x_2) - temp_max, th.tensor(0))
            C_x_3 = oc.constraint_3(x_1_scaled, x_2_scaled)  # demoulding time
            G_x_3 = c_3 * th.max(th.as_tensor(C_x_3) - time_max, th.tensor(0))
            G_x_4 = th.max(x_2 - max_agg_ratio, th.tensor(0))
            constraints = G_x_1 + G_x_2 + G_x_3 + G_x_4

            # with constraints
            c_o = 0.0001  # objective scaling
            grad_est_obj = (c_o * obj)*(q_x_1.log_prob(x_1)+q_x_2.log_prob(x_2))
            grad_est_cons = constraints*(q_x_1.log_prob(x_1)+q_x_2.log_prob(x_2))
            U_theta_holder.append(grad_est_obj+grad_est_cons)
            # w/o gradients
            #U_theta_holder.append(grad_est_obj)

            #U_theta_holder.append(
            #    (0.0001 * obj + constraints) * (q_b_1.log_prob(b_1) + q_b_2.log_prob(b_2) + q_x_2.log_prob(x_2)))

            # without constraints
            # U_theta_holder.append(obj * (q_b_1.log_prob(b_1) + q_b_2.log_prob(b_2) + q_x_2.log_prob(x_2)))
            # add to lists
            obj_holder.append(obj)
            C_1_holder.append(C_x_1)
            C_2_holder.append(C_x_2)
            C_3_holder.append(C_x_3)
        U_theta = th.sum(th.stack(U_theta_holder)) / num_samples
        with th.no_grad():
            obj_mean = np.sum(np.stack(obj_holder)) / num_samples
            C_1_mean = np.sum(np.stack(C_1_holder)) / num_samples
            C_2_mean = np.sum(np.stack(C_2_holder)) / num_samples
            C_3_mean = np.sum(np.stack(C_3_holder)) / num_samples
        assert U_theta.requires_grad == True
        return U_theta, obj_mean, C_1_mean, C_2_mean, C_3_mean

    def objective_parallel(x_1, x_2, **kwargs):
        """

        Parameters
        ----------
        x_1: aggregate,
        x_2: cem ratio
        kwargs

        Returns
        -------

        """
        if isinstance(x_2, dict):
            q_x_2 = _translate_design_variable_to_stochastic(x=x_2)

        if isinstance(x_1,dict):
            q_x_1 = _translate_design_variable_to_stochastic(x=x_1)
        num_samples = kwargs['num_samples']

        # defining holders
        U_theta_holder = []
        obj_holder = []
        C_1_holder = []
        C_2_holder = []
        C_3_holder = []

        # generate seeds

        seed_tmp = []
        X_tmp = np.ndarray(shape=(num_samples,2))
        for i in range(num_samples):
            # set seed
            # The seed will ensure that the same RV samples are passed inside the forward model
            random_seed = np.random.randint(666)
            seed_tmp.append(random_seed)
            th.manual_seed(random_seed)

            # collect RV samples
            x_2 = q_x_2.sample()
            x_1 = q_x_1.sample()

            # logistic sigmoid function to bound the input in 0-1 = 1/(1+e^(-y)), 
            # https://www.sciencedirect.com/topics/computer-science/logistic-sigmoid
            
            # TODO: ugly hardcoded, improve it
            #x_1_scaled_back = x_1.item()*(1100.0 - 700.0) + 700.0 # = x_scaled*(x_max-x_min) +x_min
            x_2_scaled = th.special.expit(x_2)

            #x_1_scaled_back = th.exp(x_1)

            x_1_scaled_back = th.special.expit(x_1).item()*(1300.0-600.0) + 600.0 # = x_scaled*(x_max-x_min) +x_min, sogmoid transform to tranform it from 0-1 and rescale to have it from 500-1300
            #X_tmp[i,0] = x_1_scaled.item()
            X_tmp[i, 0] = x_1_scaled_back  # since height need not be scaled.
            X_tmp[i,1] = x_2_scaled.item()
        # save the seed and the design varuables
        np.save('./seed_tmp.npy', np.array(seed_tmp))
        np.save('./design_var_tmp.npy',X_tmp)

        # run the workflows in parallel
        python_fn_run_jobs('../../../lebedigital/demonstrator_optimization_scripts/',no_samples=num_samples)

        for i in range(num_samples):
            th.manual_seed(seed_tmp[i])

            # collect RV samples
            x_2 = q_x_2.sample()
            x_1 = q_x_1.sample()

            # collect the kpis from the workflows
            kpi_path = '../../optimization_paper/' + str(i+1) + '/kpi.json'
            if not os.path.exists(kpi_path):
                print(f"Error: File {kpi_path} does not exist.")
                continue # some FEM solvers are not converging weirdly, so skipping those values

            obj, C_x_1, C_x_2, C_x_3 = read_kpis(kpi_path=kpi_path)

            # define constraints
            # --- Set inputs for the constraints
            # time_max = th.tensor(3)
            # temp_max = th.tensor(70)
            # max_agg_ratio = th.tensor(0.7)
            # # workability constraint. Now temp that agg ratio < 0.6
            # c_1 = 1e03
            # c_2 = 0.1
            # c_3 = 1
            # c_4 = 1
            # # design criterion
            # G_x_1 = c_1 * th.max(-th.as_tensor(C_x_1), th.tensor(0))
            # # temp
            # G_x_2 = c_2 * th.max(th.as_tensor(C_x_2) - temp_max, th.tensor(0))
            # # demoulding time
            # G_x_3 = c_3 * th.max(th.as_tensor(C_x_3) - time_max, th.tensor(0))
            # # TODO: X_tmp[i,0] below is temp for aggregate ratio.
            # #G_x_4 = th.max(th.as_tensor(X_tmp[i,0]) - max_agg_ratio, th.tensor(0))
            # constraints = G_x_1 + G_x_2 + G_x_3 #+ G_x_4

            #TODO : include the third constraint too
            #print('Including the time constraint too')
            constraints = th.max(th.as_tensor(C_x_1),th.tensor(0)) + th.max(th.as_tensor(C_x_2),th.tensor(0)) + th.max(th.as_tensor(C_x_3),th.tensor(0))

            # with constraints
            c_o = 0.001  # objective scaling (a bit less than divided by 1000. very ugly)
            grad_est_obj = (c_o * obj) * (q_x_1.log_prob(x_1) + q_x_2.log_prob(x_2))
            grad_est_cons = constraints * (q_x_1.log_prob(x_1) + q_x_2.log_prob(x_2))
            U_theta_holder.append(grad_est_obj + grad_est_cons)
            # w/o gradients
            # U_theta_holder.append(grad_est_obj)

            # U_theta_holder.append(
            #    (0.0001 * obj + constraints) * (q_b_1.log_prob(b_1) + q_b_2.log_prob(b_2) + q_x_2.log_prob(x_2)))

            # without constraints
            # U_theta_holder.append(obj * (q_b_1.log_prob(b_1) + q_b_2.log_prob(b_2) + q_x_2.log_prob(x_2)))
            # add to lists
            obj_holder.append(obj)
            C_1_holder.append(C_x_1)
            C_2_holder.append(C_x_2)
            C_3_holder.append(C_x_3)
        U_theta = th.sum(th.stack(U_theta_holder)) / num_samples
        with th.no_grad():
            U_theta_var = th.var(th.stack(U_theta_holder))
            obj_mean = np.sum(np.stack(obj_holder)) / num_samples
            C_1_mean = np.sum(np.stack(C_1_holder)) / num_samples
            C_2_mean = np.sum(np.stack(C_2_holder)) / num_samples
            C_3_mean = np.sum(np.stack(C_3_holder)) / num_samples
            # get variance of objetcive and constraints
            obj_var = np.var(np.stack(obj_holder))
            C_1_var = np.var(np.stack(C_1_holder))
            C_2_var = np.var(np.stack(C_2_holder))
            C_3_var = np.var(np.stack(C_3_holder))
        assert U_theta.requires_grad == True
        return U_theta, U_theta_var, obj_mean, C_1_mean, C_2_mean, C_3_mean, obj_var, C_1_var, C_2_var, C_3_var


    # check
    # sigma = th.tensor([1.])
    # beta = th.tensor(2 * th.log(sigma), requires_grad=True)
    # tmp = objective(x_1=th.tensor([0.4], requires_grad=True), x_2={'mean': [0.45], 'cov': th.exp(beta)}, num_samples=2)

    def optimize(design_variables:dict, eps=0.001, verbose=True, lr=0.01, number_steps: int = 10,
                 number_samples: int = 1) -> None:
        """

        Parameters
        ----------
        x1_init: slag
        x2_mu : agg ratio
        x2_sigma
        eps
        verbose
        lr
        number_steps
        number_samples

        Returns
        -------

        """
        # defining design variables
        #x_1 = th.tensor(x1_init, requires_grad=True)
        x1_mean = th.tensor(design_variables['x_1']['mean'], requires_grad=True)
        x1_sigma = th.tensor(design_variables['x_1']['s.d'], requires_grad=True)
        #beta_1 = th.tensor(2 * th.log(x1_sigma), requires_grad=True)
        x2_mean = th.tensor(design_variables['x_2']['mean'], requires_grad=True)
        x2_sigma = th.tensor(design_variables['x_2']['s.d'], requires_grad=True)
        #beta_2 = th.tensor(2 * th.log(x2_sigma), requires_grad=True)

        # defining optimizer
        parameters = [x1_mean,x1_sigma,x2_mean,x2_sigma]
        optimizer = th.optim.Adam(parameters, lr=lr)

        # value holders
        losses = []
        objective_value = []
        constraints = []
        x1_tracking = []  # Intermediate for tracking
        x2_mean_tracking = []
        x2_sigma_tracking = []
        grad_1 = []
        # Y_b_step = []
        df = pd.DataFrame()
        num_steps = number_steps
        for i in range(num_steps):
            optimizer.zero_grad()
    
            #x_1 = {'mean': x1_mean, 's.d': th.exp(0.5 * beta_1)}
            #x_2 = {'mean': x2_mean, 's.d': th.exp(0.5 * beta_2)}  # sigma = sqrt(esp(beta))
            x_1 = {'mean': x1_mean, 's.d': th.sqrt(th.exp(x1_sigma))} # passing exp transformed values to sd
            x_2 = {'mean': x2_mean, 's.d': th.sqrt(th.exp(x2_sigma))}
            loss, loss_var, obj_mean, C_1_mean, C_2_mean, C_3_mean, obj_var, C_1_var, C_2_var, C_3_var = objective_parallel(x_1=x_1, x_2=x_2, num_samples=number_samples)
            # compute grads
            loss.backward()
            # print(XX.grad)
            # losses.append(loss)
            # x1_tracking.append(x_1.clone().detach())
            # x2_mean_tracking.append(x2_mean.clone().detach())
            # # sigma_list.append(sigma)
            # x2_sigma_tracking.append(th.sqrt(th.exp(beta.clone().detach())))
            # grad_1.append(x_1.grad.clone().detach())
            if verbose:
                # if num_steps % 5 == 0:
                print(
                    f"Iteration :{i + 1}, loss value: {loss}, x1: {x1_mean}, x2_mean: {x2_mean}, sigma value: {x2_sigma}"
                    f",grad w.r.t x_1: {x1_mean.grad},,grad w.r.t x2_mean: {x2_mean.grad}"
                    f",grad w.r.t x2_sigma: {x2_sigma.grad}")


            # taking optimizer step
            optimizer.step()
            # setting bounds for the design variables
            # with th.no_grad(): # this works
            #     x1_mean.clamp_(0.1,0.8)
            #     x2_mean.clamp_(0.1,0.7) # agg ratio is set to 0.7 for workability contraints

            df = df.append({'loss': loss.item(), 'loss_var':loss_var.item(), 'objective': obj_mean, 'C_1': C_1_mean,
                              'C_2': C_2_mean, 'C_3': C_3_mean, 'x_1_mean': th.special.expit(x1_mean.clone().detach()).item(),
                            'x_1_std': np.exp(x1_sigma.clone().detach().item()),
                            #'x_1_std': np.sqrt(np.exp(beta_1.clone().detach().item())),
                            'x_2_mean': th.special.expit(x2_mean.clone().detach()).item(),
                            'x_2_std': np.exp(x2_sigma.clone().detach().item()),
                              #'x_2_std': np.sqrt(np.exp(beta_2.clone().detach().item())),
                              'x_1_mean_grad': x1_mean.grad.clone().detach().item(),
                              'x_2_mean_grad': x2_mean.grad.clone().detach().item(),
                              'obj_var': obj_var,
                              'C_1_var': C_1_var,
                              'C_2_var': C_2_var,
                              'C_3_var': C_3_var}
                             , ignore_index=True)
            df.to_csv('./Results/Optimization_results_tmp.csv',index=False)

        return df

if __name__ == '__main__':
    # x = 1/(1+e^(-y)), where y is the gaussian. so y = ln(x/(1-x)). So y mean and sd needs to be init by this.

    #x1_init = th.special.logit(th.tensor([0.25]))
    #x1_scaled_init = (1050.0 - 700.0)/(1100.0 - 700.0) # (x - x-min) / (x_max - x_min)
    # normalize the height, the range is 500-1300mm
    x1_norm = (810-600)/(1300-600)
    #x1_norm = 0.357344248553906 # if starting from some value
    x1_scaled_init = th.special.logit(th.tensor([x1_norm])) # starting from height 600 mm
    #x1_scaled_init = th.log(th.tensor([600.0])) # starting from height 600 mm
    x2_init = th.special.logit(th.tensor([0.74])) # sigmoid transformed values are passed, then later transformed back to normal.

    #x1_sigma_phi = np.log(1.0) # So that at the end we get 0.05 as exp^phi value
    #x2_sigma_phi = np.log(1.0)
    x1_sigma_phi = np.log(0.21) # So that at the end we get 0.05 as exp^phi value
    x2_sigma_phi = np.log(0.06)

    # beam height is directly proporstional to GWP, and slag ratio is inversely proportional to GWP.
    design_variables = {'x_1': {'mean': [x1_scaled_init.item()] ,'s.d': [x1_sigma_phi]},
                       'x_2': {'mean': [x2_init.item()] ,'s.d': [x2_sigma_phi]}}

    #design_variables = {'x_1': {'mean': [0.25] ,'s.d': [0.5]},
    #                    'x_2': {'mean': [0.35] ,'s.d': [0.5]}}
    df = optimize(design_variables,lr =0.1,number_steps=200,number_samples=200) # 120 step, 125 sample,
    df.to_csv('./Results/optimization_results_'+datetime+'.csv',index=False)

