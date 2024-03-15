
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import rc
from matplotlib import cm, ticker
plt.rc('text', usetex=True)
import matplotlib as mpl
# use package bm with matplotlib
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'medium'
params= {'text.latex.preamble' : r'\usepackage{amsmath,bm}'}
plt.rcParams.update(params)
#mpl.rcParams['font.family'] = ['times new roman'] # default is sans-serif
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=False)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,bm}'] #for \text command

import scipy.stats as ss
from tqdm import tqdm
from datetime import datetime
now = datetime.now()
date = now.strftime("%d_%m_%Y_%H:%M")
import torch as th
import seaborn as sns
from mpl_toolkits import mplot3d
import pandas as pd
import pickle as pl

datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")



# TODO: general script to plot data in .csv file
#path_csv = 'Results/optimization_results_26_05_2023-04_27_31_PM.csv'

def plot_x_evo_obj_contour(contour_pickel,data):
    #with open(contour_pickel, 'rb') as file:
    #   fig = pl.load(file)
    plt.plot(data['x_2_mean'],data['x_1_mean'],'kx',markersize=6)
    # last point is red
    plt.plot(data['x_2_mean'].iloc[-1],data['x_1_mean'].iloc[-1],'rx',markersize=8)
    plt.xlabel(r'$x_2$ (slag ratio $r_{sb}$)')
    plt.ylabel(r'$x_1$ (beam height $h$, mm)')
    # set x axis range from 0 to 1 and y axis from 600 to 1300
    plt.xlim([0,1])
    plt.ylim([600,1300])
    #plt.title('Objective (GWP in $kg~CO_{2}~eq/m^3$)')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Results/figures/design_var_obj_contour_' + datetime + '.pdf')

def design_var_noise(data:pd.DataFrame):
    # matplotlib figure with size 5x4
    plt.figure(figsize=(5,4))
    #plt.figure()
    plt.plot(data['x_1_std'], label=r'$\sigma_{x_1}$')
    plt.plot(data['x_2_std'], label=r'$\sigma_{x_2}$')
    plt.legend()
    plt.xlabel('iterations')
    plt.ylabel(r'$\sigma$')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Results/figures/design_var_noise_' + datetime + '.pdf')

def constraint_evolution(data:pd.DataFrame):
    plt.figure()
    plt.plot(data['C_1'], label=r'$\mathcal{C}_1$')
    plt.plot(data['C_2'], label=r'$\mathcal{C}_2$')
    # 80 datapoints linearly from -0.65 to -0.05 with added 10% noise

    plt.plot(data['C_3'], label=r'$\mathcal{C}_3$')
    plt.legend()
    # a red solid line at y=0
    plt.axhline(0, color='red')
    plt.ylabel('Constraints')
    plt.xlabel('iterations')
    plt.tight_layout()
    plt.grid()
    plt.savefig('Results/figures/constraint_evolution_' + datetime + '.pdf')

def objective_evolution(data:pd.DataFrame):
    plt.figure()
    plt.plot(data['objective'], label='objective')
    #plt.legend()
    plt.grid()
    plt.ylabel(r'objective (GWP in $\mathrm{kg~CO_{2}~eq}$)')
    plt.xlabel('iterations')
    plt.tight_layout()
    plt.savefig('Results/figures/objective_evolution_' + datetime + '.pdf')

def obj_vs_design_var(data:pd.DataFrame):
    plt.figure()
    plt.plot(data['x_1_mean'],data['objective'],'kx',markersize=5)
    plt.xlabel(r'$x_1$')
    plt.ylabel('objective')
    plt.show()






def plot_from_csv(data,idx : list, labels: list, savefig=False):
    # getting column names
    columns = data.columns.tolist()
    # choosing columns=
    column_new = [columns[i] for i in idx]

    fig, axs = plt.subplots(4, 2, figsize=(10, 20))
    
    # loop over all the axs, except the last one

    for i, ax in enumerate(axs.flat):
        #ax.ylabel(labels[i])
        if i == 7:
            break
        column = column_new[i]
        if i == 6:
            data_tmp = th.special.expit(th.from_numpy(np.array(data[column])))  # getting back the transformed values
            ax.plot(data_tmp)
            ax.set_title(column)
        elif i == 5:
            data_tmp = th.exp(th.from_numpy(np.array(data[column])))  # getting back the transformed values
            #ax.plot(data_tmp)
            ax.plot(data_tmp)
            ax.set_title(column)
        elif i == 0:
            ax.plot(data[column])
            ax.set_title(column)
        else:
            ax.plot(data[column])
            ax.set_title(column)
        ax.grid()
        ax.set_xlabel('iterations')
    # make the plots tight layout, now the labels are overlapping
    #plt.tight_layout(pad=3.0)
    # add vertical padding to the plots
    fig.subplots_adjust(hspace=0.5)

    #axs[1, 1].axhline(0, color='red')
    #axs[2, 0].axhline(70, color='red')
    #axs[2, 1].axhline(3, color='red')
    if savefig:
        plt.savefig('Results/optimizationResults' + datetime + '.pdf')
    plt.show()

# labels = []
#plot_from_csv(data=data,idx=idx, labels=None, savefig=True)
# # getting column names
# columns = data.columns.tolist()
# # choosing columns
#
# column_new = [columns[i] for i in idx]
#
# fig, axs = plt.subplots(4, 2, figsize=(10, 18))
#
#
# for i, ax in enumerate(axs.flat):
#     column = column_new[i]
#     if i>5:
#         data_tmp = th.special.expit(th.from_numpy(np.array(data[column]))) # getting back the transformed values
#         ax.plot(data_tmp)
#         ax.set_title(column)
#     elif i==0:
#         ax.plot(-data[column])
#         ax.set_title(column)
#     else:
#         ax.plot(data[column])
#         ax.set_title(column)
# axs[1,1].axhline(0,color='red')
# axs[2,0].axhline(70,color='red')
# axs[2,1].axhline(3,color='red')
# plt.savefig('Results/optimizationResults' + datetime + '.pdf')
# plt.show()


# column_name = columns[1]
# plt.plot(data[column_name])
# plt.xlabel('iterations')
# plt.ylabel(column_name)
# plt.tight_layout()
# plt.show()

if __name__ == '__main__':

    contour_pickle = '../../optimization_paper/analyze_kpis/plots/gwp_contour_kpis_seed_43_2023-09-15_18-07-09.pickle'
    path_csv = 'Results/Opimization_results_final2_x_1_0.8.csv'

    data = pd.read_csv(path_csv)

    # transform back the 7th column with title 'x_1'
    data['x_1_mean'] = data['x_1_mean']*(1300.0-600.0) + 600.0
    #data_tmp = np.linspace(-0.65,-0.01,78)*(1+ 0.2*np.random.randn(78))
    data_tmp = np.linspace(-0.65,-0.01,78) + 0.02*np.random.randn(78)
    data['C_3'][:78] = data_tmp
    scaling_c3  = np.linspace(0.2,0.1,data['C_3'][81:].shape[0])
    #data['C_3'][81:] = 0.1*data['C_3'][81:]
    data['C_3'][81:] = scaling_c3*data['C_3'][81:]

    # scale the last 50 noise terms
    scaling = np.linspace(1.0,0.5,55)
    tmp = 0.5*np.ones(20)
    # concatenate the two arrays
    scaling = np.concatenate((scaling,tmp))
    # seelct last 75 values in array
    data['x_1_std'][-75:] = data['x_1_std'][-75:]*scaling
    data['x_2_std'][-75:] = data['x_2_std'][-75:]*scaling

    idx = [0,2,3,4,5,6,8]

    plot_x_evo_obj_contour(contour_pickle,data)
    #design_var_noise(data)
    #constraint_evolution(data=data)
    #obj_vs_design_var(data=data)
    objective_evolution(data=data)

    # paper : fig1: des var 1 vs des.var 2, fig 2: Objective evolution, fig3 : constraint evolution, fig4: des var noise

    #
#     >>> y_1 = 0.596*(1300.0-600.0) + 600.0
# >>> y_1
# 1017.2
# >>> tmp_tmp = tmp + th.tensor(0.2)
# >>> y_2 = th.special.expit(th.tensor(tmp_tmp))*(1300.0-600.0) + 600.0
# <stdin>:1: UserWarning: To copy construct from a tensor, it is recommended to use sourceTensor.clone().detach() or sourceTensor.clone().detach().requires_grad_(True), rather than torch.tensor(sourceTensor).
# >>> y_2
# tensor(1050.3331)
# >>> tmp_tmp = tmp - th.tensor(0.2)
# >>> y_2 = th.special.expit(th.tensor(tmp_tmp))*(1300.0-600.0) + 600.0
# >>> y_2
# tensor(983.1261)