import torch as th
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import copy
import pandas as pd
import csv
import sys, pathlib

from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
import matplotlib as mpl
from matplotlib import rc

from lebedigital.demonstrator_calibration.parametric_model import NN_mean, train_NN
from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.demonstrator_calibration.likelihood import gaussian_likelihood
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data
from lebedigital.demonstrator_calibration.sampler import MCMC_DRAM
from lebedigital.demonstrator_calibration.VBEM import VBEM


# 25.Aug 4pm, restarted as the simulation crashed after 250 iterations. restarted from that point.  

def main():
    hydration_solver = HydrationSolverWrapper()
    #cov_param = th.tensor([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],requires_grad=True)
    cov_param = th.tensor([0.5,0.0,0.5,0.0,0.0,0.5,0.0,0.0,0.0,0.5],requires_grad=True)
    #cov_param = th.tensor([1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],requires_grad=True)
    #data_path = 'usecases/optimization_paper/calibration_data/Excel_files/hydration_data_processed.xlsx'
    data_path = '../../../calibration_data/Excel_files/hydration_data_processed.xlsx'
    df_data = process_hydration_data(data_path)
    b_opt = np.load('../../../calibration_data/optimization_results_hydration.npy')
    #b_opt = np.load('usecases/optimization_paper/calibration_data/optimization_results_hydration.npy')

    # log transfor the B2 as scaling
    b_opt[:,1] = b_opt[:,1]*1e-03
    b_opt[:,1] = np.log(b_opt[:,1])
    b_init = b_opt.tolist()

    # add noise to the b_init
    #b_opt_log_ = th.tensor(b_opt_log)*(1 + th.randn(b_opt_log.shape)*0.01) # adding 10% of the data as noise
    # mean = np.mean(b_opt,axis=0)
    # std = np.std(b_opt,axis=0)
    # b_opt = (b_opt-mean)/std
    b_opt = th.tensor(b_opt)*(1 + th.randn(b_opt.shape)*0.1) # adding noise
    #print(f'b_opt = {b_opt}')
    #breakpoint()

    # the simulation crashed somehow, so reastarting from that point
    restart_from_a_point = True
    if restart_from_a_point:
        cov_param = np.genfromtxt('cov_parameters2023_08_24-03_54_27_PM.csv', delimiter=',').tolist()[-1]
        cov_param = th.tensor(cov_param,requires_grad=True)
        nn_mean = NN_mean(input_dim=1, output_dim=4, hidden_dim=20)
        nn_mean.load_state_dict(th.load('NN_state_dict_till_itr_200_2023_08_24-03_54_27_PM.pth'))
        b_init = [[3.0480967, -8.4496829,  3.238852 ,  6.9784307],[4.8550468, -13.793785 ,   4.5960467,   4.0132815]
                  ,[3.522844 , -8.0200439,  5.5433366,  2.5466527],[0.08126229, -0.31569943,  1.0669216 ,  1.130067]]

        vbem = VBEM(prior=prior,forward_model=hydration_solver.solve,likelihood=gaussian_likelihood,
                    model_prior_mean=nn_mean,prior_cov_params=cov_param, sigma_likelihood=3, latent_dim=4,
                    dataframe_observed_data=df_data,no_observed_data_pair=4,b_init=b_init,
                    pre_train=False,lr=5e-03) # smaller lr to avoid divergence
    else:    
        vbem = VBEM(prior=prior,forward_model=hydration_solver.solve,likelihood=gaussian_likelihood,
                    model_prior_mean=NN_mean,prior_cov_params=cov_param, sigma_likelihood=3, latent_dim=4,
                    dataframe_observed_data=df_data,no_observed_data_pair=4,b_init=b_init,
                    pre_train=True,lr=1e-02)
    
    vbem.M_step(201,no_samples=100)

if __name__ == '__main__':
    # add the path of the current file
    sys.path.append(str(pathlib.Path(__file__).resolve().parent))
    # assert that this function is called from its parent
    main()

