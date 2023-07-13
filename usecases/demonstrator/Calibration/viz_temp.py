
import numpy as np
import matplotlib.pyplot as plt
plt.style.use({'figure.facecolor':'white'})
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib import rc
from matplotlib import cm, ticker
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

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

datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")



# TODO: general script to plot data in .csv file
path_csv = 'Results/optimization_results_26_05_2023-04_27_31_PM.csv'

data = pd.read_csv(path_csv)
idx = [0,1,2,3,4,5,6,7,8,9]

def plot_from_csv(data,idx : list, labels: list, savefig=False):
    # getting column names
    columns = data.columns.tolist()
    # choosing columns=
    column_new = [columns[i] for i in idx]

    fig, axs = plt.subplots(5, 2, figsize=(10, 22))

    for i, ax in enumerate(axs.flat):
        #ax.ylabel(labels[i])
        column = column_new[i]
        if i == 6 or i == 8:
            data_tmp = th.special.expit(th.from_numpy(np.array(data[column])))  # getting back the transformed values
            ax.plot(data_tmp)
            ax.set_title(column)
        elif i == 0:
            ax.plot(-data[column])
            ax.set_title(column)
        else:
            ax.plot(data[column])
            ax.set_title(column)
    axs[1, 1].axhline(0, color='red')
    axs[2, 0].axhline(70, color='red')
    axs[2, 1].axhline(3, color='red')
    if savefig:
        plt.savefig('Results/optimizationResults' + datetime + '.pdf')
    plt.show()

# labels = []
plot_from_csv(data=data,idx=idx, labels=None, savefig=True)
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

print(i)
# column_name = columns[1]
# plt.plot(data[column_name])
# plt.xlabel('iterations')
# plt.ylabel(column_name)
# plt.tight_layout()
# plt.show()