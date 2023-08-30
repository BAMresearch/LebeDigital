import numpy as np
import torch as th
from matplotlib import pyplot as plt
import seaborn as sb
# use latex with matplotlib
plt.rc('text', usetex=True)
import matplotlib as mpl
params= {'text.latex.preamble' : r'\usepackage{amsmath,bm}'}
plt.rcParams.update(params)
# mpl.rcParams['font.size'] = 14
# mpl.rcParams['legend.fontsize'] = 'medium'

from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.demonstrator_calibration.parametric_model import NN_mean
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data, process_homogenization_data
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper, HomogenizationSolverWrapper
from lebedigital.demonstrator_calibration.visualization import plot_data, viz_learnt_prior_model, prob_hydration_solver_output

import sys, pathlib
from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")



# plot learnt prior model
nn_model = NN_mean(input_dim=1, output_dim=4, hidden_dim=20)
nn_state_dict = 'usecases/optimization_paper/model_learning/hydration/exp_11/NN_state_dict_till_itr_200_2023_08_25-04_41_31_PM.pth'
cov_path = 'usecases/optimization_paper/model_learning/hydration/exp_11/cov_parameters2023_08_24-03_54_27_PM.csv'
cov_params = np.genfromtxt(cov_path, delimiter=',').tolist()[-1]
def transformed_back(samples):
    shape = samples.shape
    # exp transform to the last dimention
    samples[:,:, 0] = samples[:,:, 0] * 1e-04
    samples[:,:, 1] = np.exp(samples[:,:, 1])
    samples[:,:, 3] = samples[:,:, 3] * 1e05
    assert samples.shape == shape, "shape of the samples is changed"
    return samples

viz_learnt_prior_model(nn_model,nn_state_dict,cov_params,latent_dim=4,case='hydration',transform_unscaled=transformed_back,
                       save_path='usecases/optimization_paper/model_learning/hydration/exp_11/Results/')#,save_path='lebedigital/demonstrator_calibration/misc/')

# clip values of the list to 0.1 if its more than that
# for i in range(len(cov_params)):
#     if cov_params[i] > 0.1:
#         cov_params[i] = 0.1

# plot the evolution of the loss


# plot evolution of cov parameters
legend = [r'$\phi_{11}$',r'$\phi_{21}$',r'$\phi_{22}$',r'$\phi_{31}$',r'$\phi_{32}$',r'$\phi_{33}$',r'$\phi_{41}$',r'$\phi_{42}$',r'$\phi_{43}$',r'$\phi_{44}$']
plot_data(path=cov_path, labels=[r'$\bm{\phi$}'],legends=legend,save_path='usecases/optimization_paper/model_learning/hydration/exp_11/Results/cov_parameters'+datetime+'.pdf')

prob_hydration_solver_output(NN_model=nn_model, NN_state_dict=nn_state_dict, cov_params=cov_params, latent_dim=4, 
                             save_path='usecases/optimization_paper/model_learning/hydration/exp_11/Results/')

