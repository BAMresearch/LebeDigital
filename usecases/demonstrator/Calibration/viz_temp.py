
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
path_csv = 'Results/optimization_results_25_04_2023-01_42_05_PM.csv'

data = pd.read_csv(path_csv)

# getting column names
columns = data.columns.tolist()
# choosing columns
idx = [1,2,3,4,5,7]
column_new = [columns[i] for i in idx]

fig, axs = plt.subplots(3, 2, figsize=(10, 12))


for i, ax in enumerate(axs.flat):
    column = column_new[i]
    if i>3:
        data_tmp = th.special.expit(th.from_numpy(np.array(data[column]))) # getting back the transformed values
        ax.plot(data_tmp)
        ax.set_title(column)
    else:
        ax.plot(data[column])
        ax.set_title(column)
axs[0,1].axhline(0,color='red')
axs[1,0].axhline(70,color='red')
axs[1,1].axhline(3,color='red')
plt.savefig('Results/optimizationResults' + datetime + '.pdf')
plt.show()

print(i)
# column_name = columns[1]
# plt.plot(data[column_name])
# plt.xlabel('iterations')
# plt.ylabel(column_name)
# plt.tight_layout()
# plt.show()