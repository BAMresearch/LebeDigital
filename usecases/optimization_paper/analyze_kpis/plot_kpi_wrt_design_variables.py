#%%
import numpy as np
import torch as th
from matplotlib import pyplot as plt
import seaborn as sb
import pandas as pd
# use latex with matplotlib
plt.rc('text', usetex=True)
import matplotlib as mpl
# use package bm with matplotlib
mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 'medium'
params= {'text.latex.preamble' : r'\usepackage{amsmath,bm}'}
plt.rcParams.update(params)

import time
datetime = time.strftime("%Y-%m-%d_%H-%M-%S")

#%%
def kpi_vs_x(csv_file:str, combined:bool =True):
    # load the csv file
    data = pd.read_csv(csv_file)
    fig, ax = plt.subplots(1,1)
    for col in data.columns[2:]:
        if combined:
            dim = np.unique(data['height'].values).shape[0]
            if col == 'gwp':
                gwp = ax.contourf(np.unique(data['slag_ratio'].values),np.unique(data['height'].values),data[col].values.reshape(dim,dim))
                fig.colorbar(gwp)
            elif col == 'constraint_beam_design':
                # plot of single sontour line as an indicator fucntion
                ax.contour(np.unique(data['slag_ratio'].values),np.unique(data['height'].values),data[col].values.reshape(dim,dim),levels=[0],colors='k', linestyles='dashed')
            elif col == 'constraint_temperature':
                ax.contour(np.unique(data['slag_ratio'].values),np.unique(data['height'].values),data[col].values.reshape(dim,dim),levels=[0],colors='r', linestyles='dashed')
            elif col == 'constraint_time':
                ax.contour(np.unique(data['slag_ratio'].values),np.unique(data['height'].values),data[col].values.reshape(dim,dim),levels=[0],colors='b', linestyles='dashed')
            ax.set_xlabel(r'$x_2$')
            ax.set_ylabel(r'$x_1$')
            # title of the plot

            ax.set_title('GWP and constraints')
            # save the plot
            plt.savefig(f'plots/combined_contour.pdf')
        # plot the contours
        else:
            plt.figure()
            dim = np.unique(data['height'].values).shape[0]
            plt.contourf(np.unique(data['slag_ratio'].values),np.unique(data['height'].values),data[col].values.reshape(dim,dim))
            # colorbar for the above
            plt.colorbar()
            plt.title(f'{col}')
            # set colorbar for the axis
            #plt.colorbar()
            plt.xlabel(r'$x_2$')
            plt.ylabel(r'$x_1$')
            plt.savefig(f'plots/{col}_contour.pdf')



    

#%%
# Update 1st Sep,2023. Found a set of inputs which gives a good optimization problem. See the latest plot. Obj/constraint variability needs to be checked.
csv_file = 'kpis_2023-08-31_17-28-38.csv'

kpi_vs_x(csv_file,combined=True)

#%%
