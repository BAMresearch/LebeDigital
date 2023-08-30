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
from lebedigital.demonstrator_calibration.forward_solvers import HomogenizationSolverWrapper
from usecases.optimization_paper.calibration_data.data_handling import process_homogenization_data
from lebedigital.demonstrator_calibration.sampler import MCMC_DRAM
from lebedigital.demonstrator_calibration.VBEM_homogenization import VBEM

# this exp with pre train with 1000 steps

def main():
    homogenization_solver = HomogenizationSolverWrapper()
    cov_param = th.tensor([1.0,0.0,1.0],requires_grad=True)
    #data_path = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'
    data_path = '../../../calibration_data/Excel_files/homogenization_data_processed_E.csv'
    df_data = process_homogenization_data(data_path)
    #b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_2023_08_21-05_47_52_PM.npy')
    b_opt = np.load('../../../calibration_data/optimization_results_homogenization_2023_08_21-05_47_52_PM.npy')
    #b_opt = np.delete(b_opt,3,0) # deleting the 3rd row
    b_opt[:,0] = b_opt[:,0]*1e-09
    b_opt[:,1] = b_opt[:,1]*1e-07
    b_opt = th.tensor(b_opt)*(1 + th.randn(b_opt.shape)*0.05) # adding noise
    # create new tensor by removing the 4th row
    #b_opt = th.cat((b_opt[:,0:3],b_opt[:,4:]),dim=1)

    b_init = b_opt.tolist()
    restart_from_a_point = True
    if restart_from_a_point:
        cov_param = np.genfromtxt('cov_parameters2023_08_27-03_04_50_PM.csv', delimiter=',').tolist()[-1]
        cov_param = th.tensor(cov_param,requires_grad=True)
        nn_mean = NN_mean(input_dim=1, output_dim=2, hidden_dim=20)
        nn_mean.load_state_dict(th.load('NN_state_dict_till_itr_150_2023_08_27-03_04_50_PM.pth'))
        b_init = [[9.7285462,  3.2614563],[9.5320265,  3.2381915],[9.7427698, 3.1489173],
                  [8.6139817, 2.6081752],[8.1403589, 2.2345528 ],[7.2903925, 1.7900311]]

        vbem = VBEM(prior=prior,forward_model=homogenization_solver.solve,likelihood=gaussian_likelihood,
                model_prior_mean=nn_mean,prior_cov_params=cov_param, sigma_likelihood=[2e09,2e06], latent_dim=2,
                dataframe_observed_data=df_data,no_observed_data_pair=6,b_init=b_init,
                pre_train=False,lr=5e-03)
    else:
        vbem = VBEM(prior=prior,forward_model=homogenization_solver.solve,likelihood=gaussian_likelihood,
                model_prior_mean=NN_mean,prior_cov_params=cov_param, sigma_likelihood=[2e09,2e06], latent_dim=2,
                dataframe_observed_data=df_data,no_observed_data_pair=6,b_init=b_init,
                pre_train=True,lr=1e-02)
    
    #nn_model = NN_mean(input_dim=1, output_dim=4, hidden_dim=20)
    vbem.M_step(201,no_samples=100)

if __name__ == '__main__':
    # add the path of the current file
    sys.path.append(str(pathlib.Path(__file__).resolve().parent))
    # assert that this function is called from its parent
    main()

