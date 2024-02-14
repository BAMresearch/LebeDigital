#%%
import numpy as np
import torch as th
import matplotlib.pyplot as plt
import os
import yaml
import seaborn as sb
import pandas as pd
from datetime import datetime
import sys, pathlib
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
plt.rc('text', usetex=True)
import matplotlib as mpl
# use package bm with matplotlib
mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 'medium'
# add preamble to use latex in matplotlib
params= {'text.latex.preamble' : r'\usepackage{amsmath,bm}'}
plt.rcParams.update(params)

# local imports
from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper, HomogenizationSolverWrapper
from lebedigital.demonstrator_scripts.youngs_modulus_approximation import youngs_modulus_approximation
# import fenics_concrete

# pandas to read from excel /
# it should read the first three 3 row headers
#%%
file_location = os.path.dirname(os.path.realpath(__file__))
path = file_location + '/Excel_files/hydration_data_processed.xlsx'
#df = pd.read_excel(path, header=[0,1,2])

# path_yaml = 'usecases/optimization_paper/calibration_data/hydration_data_KIT.yaml'
# def read_yml(path):
#     with open(path) as file:
#         hydration_data = yaml.safe_load(file)
#     return hydration_data

# hy_data = read_yml(path_yaml)

def process_hydration_data(path:str, array:bool=False):
    """read the excel data, save in dict, convert time units from hours to sec and 
    convert heat to a cumulative sum

    Parameters
    ----------
    path : str
        relative path of the hydration data in .xlxs. Temp fix

    Returns
    -------
    df : dict
        dictionary of the hydration data
    """
    df = pd.read_excel(path, header=[0,1,2])

    # convert to dictionary without the index in pandas dataframe
    df = df.to_dict(orient='list')

    # drop NAN from the dictionary
    for key in df.keys():
        df[key] = [x for x in df[key] if str(x) != 'nan']

    # the key 'Age' should be converted to seconds from hours
    for k in df.keys():
        if 'Age' in k:
            df[k] = np.array([x*3600 for x in df[k]])

    # the values in key 'Q' should be converted to cumulative sum
    for k in df.keys():
        if 'Q' in k:
            df[k] = np.cumsum(df[k])
    return df

def process_homogenization_data(path_to_csv, generate_E:bool=False):
    """ output dict with keys 'x', 'fc(Mpa)' and 'E' with values as lists"""
    # read the csv file with pandas
    df = pd.read_csv(path_to_csv, header=[0])
    # convert to dictionary without the index in pandas dataframe
    df = df.to_dict(orient='list')
    if generate_E:
        density = 2820 # kg/m3, value obtained by running the homogenization code

        E = youngs_modulus_approximation(fc = df['fc(Mpa)']*ureg('MPa'), density = density*ureg('kg/m^3')).magnitude

        # add E to the dict df
        df['E'] = E.tolist() # E in MPa too

        # save the dict as a csv file
        df = pd.DataFrame.from_dict(df)
        path_to_csv = path_to_csv[:-4] + '_E.csv'
        df.to_csv(path_to_csv, index=False)

    return df


def df_to_dict_hydration(df:pd.DataFrame):
    observed = {}
    observed['x'] = np.array([[0.0],[0.3],[0.5],[0.85]])
    # convert the array to list
    observed['z'] = np.array([df[('20C','CP0','Q')],df[('20C','CP30','Q')],df[('20C','CP50','Q')],df[('20C','CP85','Q')]])
    # stack columns of dataframe to forma 2d array


    return observed

def input_hydration_solver(df:pd.DataFrame)->dict:
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = [0,5000,10000,20000,100000]
    return inp_solver


def plot_hydration_data(df):
    """To plot the hydration data for the diff ratios

    Parameters
    ----------
    df : pandas dataframe
        hydration data
    """
    # plot the data

    if len(df) == 24: # it various temp are considered
        fig, axs = plt.subplots(1,3, figsize=(13.5,4.5))
        fig.tight_layout(pad=2.0)
        ratio_keys = ['CP0','CP30','CP50','CP85']
        temp_keys = ['10C','20C','35C']
        for i in range(len(temp_keys)):
            axs[i].set_title(temp_keys[i])
            for j in range(len(ratio_keys)):
                axs[i].plot(df[(temp_keys[i],ratio_keys[j],'Age')], df[(temp_keys[i],ratio_keys[j],'Q')],'+', 
                            label=ratio_keys[j])
            #axs[i].legend()
        axs[2].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
        #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            # the legend goes outside the plot, so we need to make the plot a bit wider
        plt.subplots_adjust(right=0.90)
    else:
        fig, ax = plt.subplots(1,1)

        ratio_keys = ['CP0','CP30','CP50','CP85']
        for i in range(len(ratio_keys)):
            ax.plot(df[('20C',ratio_keys[i],'Age')], df[('20C',ratio_keys[i],'Q')],'+', 
                    label=ratio_keys[i])
        ax.legend()
        #ax.set_xlabel('Age (s)')
        #ax.set_ylabel('Cum. Heat of hydration (J/gh)')
    for ax in axs.flat:
        ax.grid()
        ax.set(xlabel='Age (s)', ylabel=r'Cum. Heat of hydration $\bm{Q}$ (J/gh)')
        ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    #fig.savefig('usecases/optimization_paper/calibration_data/figs/hydration_data' + datetime + '.png')
    fig.savefig('figs/hydration_data' + datetime + '.pdf')
    return fig

def compare_hydration_data_solver(df):
        # -- observed inputs
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = df[('20C','CP0','Age')]

    # -- latents -----
    b = [2.916,2.4229,5.554,5]
    hydration_solver_wrapper = HydrationSolverWrapper()
    heat_list = hydration_solver_wrapper(b,inp_solver)

    fig, ax = plt.subplots(1,1)

    ratio_keys = ['CP0','CP30','CP50','CP85']
    for i in range(len(ratio_keys)):
        ax.plot(df[('20C',ratio_keys[i],'Age')], df[('20C',ratio_keys[i],'Q')],'*-', label=ratio_keys[i]+'exp')
    ax.plot(df[('20C','CP0','Age')], heat_list,'X-', label='CP0sim')
    ax.legend()
    ax.set_xlabel('Age (s)')
    ax.set_ylabel('Cum. Heat of hydration (J/gh)')
    fig.savefig('usecases/optimization_paper/calibration_data/hydration_data_solver_comparision' + datetime + '.png')
    return fig

if __name__ == '__main__':

    #df = process_hydration_data(path)

    file_location = os.path.dirname(os.path.realpath(__file__))
    path_to_csv = file_location + '/Excel_files/homogenization_data_processed.csv'
    #dict_hydration = process_homogenization_data(path_to_csv=path_to_csv,generate_E=True)

    # function to plot the hydration data, age vs heat of hydration for different w/c ratios
    #%%
    
    


    # fig = plot_hydration_data(df)

    # fig_2 = compare_hydration_data_solver(df)


    # print(df)
    path_hydration_data = file_location + '/Excel_files/hydration_data_processed.xlsx'
    df = process_hydration_data(path_hydration_data)

    plot_hydration_data(df)

    #obs = df_to_dict_hydration(df)
    #print(obs)


    




# %%
