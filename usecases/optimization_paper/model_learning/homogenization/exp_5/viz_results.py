import numpy as np
import torch as th
from matplotlib import pyplot as plt
import seaborn as sb
# use latex with matplotlib
# plt.rc('text', usetex=True)
# import matplotlib as mpl
# # use package bm with matplotlib
# mpl.rcParams['font.size'] = 14
# mpl.rcParams['legend.fontsize'] = 'medium'
# params= {'text.latex.preamble' : r'\usepackage{amsmath,bm}'}
# plt.rcParams.update(params)

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman']})
import matplotlib as mpl    
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# plt.style.use('ggplot')
SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 20

rc('font', size=MEDIUM_SIZE)          # controls default text sizes
rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure titles

from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.demonstrator_calibration.parametric_model import NN_mean
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data, process_homogenization_data
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper, HomogenizationSolverWrapper
from lebedigital.demonstrator_calibration.visualization import plot_data, viz_learnt_prior_model, prob_homogenization_solver_output

import sys, pathlib
from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

path_to_EM_results = 'usecases/optimization_paper/model_learning/homogenization/exp_5/EM_results2023_08_27-03_04_50_PM.csv'
path_to_cov = 'usecases/optimization_paper/model_learning/homogenization/exp_5/cov_parameters2023_08_27-03_04_50_PM.csv'
data = np.genfromtxt(path_to_EM_results,delimiter=',')
data_cov = np.genfromtxt(path_to_cov,delimiter=',')

# plt.figure()
# # tight layouut
# plt.tight_layout()
# plt.plot(data[:,0])
# plt.xlabel('Iteration')
# plt.ylabel('$loss$')
# plt.savefig('usecases/optimization_paper/model_learning/homogenization/exp_5/Results/EM_results.pdf')
# plt.show()

def transformed_back(samples):
    shape = samples.shape
    # exp transform to the last dimention
    samples[:,:, 0] = samples[:,:, 0] * 1e09
    samples[:,:, 1] = samples[:,:, 1]* 1e07
    assert samples.shape == shape, "shape of the samples is changed"
    return samples

legends = [r'$\phi_{11}$',r'$\phi_{21}$',r'$\phi_{22}$']
#plot_data(path=path_to_cov,labels=[r'$\bm{\phi$}'],legends=legends,save_path='usecases/optimization_paper/model_learning/homogenization/exp_5/Results/cov_parameters'+datetime+'.pdf')

nn_model = NN_mean(input_dim=1, output_dim=2, hidden_dim=20)
nn_state_dict = 'usecases/optimization_paper/model_learning/homogenization/exp_5/NN_state_dict_till_itr_150_2023_08_27-03_04_50_PM.pth'
#cov_path = 'usecases/optimization_paper/model_learning/hydration/exp_11/cov_parameters2023_08_24-03_54_27_PM.csv'
cov_params = np.genfromtxt(path_to_cov, delimiter=',').tolist()[-1]
viz_learnt_prior_model(nn_model,nn_state_dict,cov_params,latent_dim=2,transform_unscaled=transformed_back,
                        case='homogenization',save_path='usecases/optimization_paper/model_learning/homogenization/exp_5/Results/')#,save_path='lebedigital/demonstrator_calibration/misc/')

prob_homogenization_solver_output(nn_model,nn_state_dict,cov_params,latent_dim=2,save_path='usecases/optimization_paper/model_learning/homogenization/exp_5/Results/')