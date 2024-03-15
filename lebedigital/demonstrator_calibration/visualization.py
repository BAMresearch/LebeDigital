import numpy as np
import torch as th
from matplotlib import pyplot as plt
import seaborn as sb
# use latex with matplotlib
plt.rc('text', usetex=True)
import matplotlib as mpl
# use package bm with matplotlib
mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 'medium'
import matplotlib as mpl
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.demonstrator_calibration.parametric_model import NN_mean
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data, process_homogenization_data
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper, HomogenizationSolverWrapper

import sys, pathlib
from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

# function to take a path to a csv file, adn the keys for the columns of interest, and plot the data
def plot_data(path:str, labels:list=None,legends:list=None,save_path=None):
    # open csv file and read data as numpy array
    data = np.genfromtxt(path, delimiter=',')
    # select only 150 rows
    #fig, ax = plt.subplots(1,1)
    # plot all the columns
    plt.figure()
    if len(labels) == 1:
        for i in range(data.shape[1]):
            # the x axis of the plot should be the index of the data
            if legends is not None:
                plt.plot(np.linspace(0,data.shape[0],data.shape[0]),data[:,i], label=legends[i])
            else:
                plt.plot(np.linspace(0,data.shape[0],data.shape[0]),data[:,i])
        # set labels
        plt.xlabel('$iterations$')
        plt.ylabel(labels[0])

        # show figure
        #plt.show()
        # set legend
        #plt.legend()
        # legend outside the plot, to the right side
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # the legend goes outside the plot, so we need to make the plot a bit wider
        plt.subplots_adjust(right=0.8)

        # save figure
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()

    # create syuplot for each column corresponding to each label
    else:
        # create subplots with no of rows and columns depending on the number of labels
        num_columns = data.shape[1]
        num_rows = (num_columns + 1) // 2  # Calculate number of subplot rows

        fig, axs = plt.subplots(num_rows, 2, figsize=(12, 6))
        
        # make the plots tight
        fig.tight_layout(pad=3.0)
        # iterate over the labels
        for i in range(len(labels)):
            row = i // 2
            col = i % 2
            # plot the column corresponding to the label
            axs[row,col].plot(np.linspace(0,data.shape[0],data.shape[0]),data[:,i])
            # set x label
            axs[row,col].set_xlabel('$iterations$')
            # set y label
            axs[row,col].set_ylabel(labels[i])
        # If there's an odd number of columns, hide the last subplot
        if num_columns % 2 == 1:
            fig.delaxes(axs[num_rows - 1, 1])
            # centre the last subplot which is visible
            #axs[num_rows - 1, 0].set_position([0.125, 0.1, 0.775, 0.8])
        # save figure
        if save_path is not None:
            plt.savefig(save_path)
        # show figure
        plt.show()

def viz_learnt_prior_model(NN_model:object,NN_state_dict:str,cov_params:list,latent_dim:int,case:str,transform_unscaled:callable = None,
                           save_path=None):
    """_summary_

    Parameters
    ----------
    NN_model : object
        _description_
    NN_state_dict : str
        _description_
    cov_pararms : list
        np.genfromtxt(path to csv file, delimiter=',').tolist(). this should be the last value/converged value
    latent_dim : int
        _description_
    case : str
        "hydration" or "homogenization"
    transform_unscaled : callable, optional
        pass a function to scale back
    save_path : _type_, optional
        _description_, by default None
    """
    # load the state dictionary
    NN_model.load_state_dict(th.load(NN_state_dict))

    # covert the values in csv file to a list
    #cov_params = np.genfromtxt(cov_pararms, delimiter=',').tolist()

    prior_model = prior(mean=NN_model, cov_params=cov_params
                           , cov_type='full',latent_dim=latent_dim) # pass the last value of the cov_params
    x_test = th.arange(0,1.0,0.01).reshape(-1,1)
    pred_sample =prior_model.sample(x_test,n_samples=500)
    if transform_unscaled is not None:
        pred_sample = transform_unscaled(pred_sample)
    # mean along the rows or the 0th dimention
    b_mean = np.mean(pred_sample,axis=0)
    b_std = np.std(pred_sample,axis=0)
    x_test = x_test.detach().numpy()

    if case == 'hydration':
        fig, axs = plt.subplots(2, 2)
            # make the plots tight
        fig.tight_layout(pad=2.0)
        axs[0, 0].plot(x_test, b_mean[:,0])
        axs[0,0].fill_between(x_test.ravel(), b_mean[:,0] - 3*b_std[:,0], b_mean[:,0] + 3*b_std[:,0], alpha=0.3)
        axs[0, 0].set_ylabel(r'$B_1, \mathrm{1/s}$')
        axs[0, 1].semilogy(x_test, b_mean[:,1])
        axs[0,1].fill_between(x_test.ravel(), b_mean[:,1] - 3*b_std[:,1], b_mean[:,1] + 3*b_std[:,1], alpha=0.3)
        axs[0, 1].set_ylabel('$B_2$')
        axs[1, 0].plot(x_test, b_mean[:,2])
        axs[1,0].fill_between(x_test.ravel(), b_mean[:,2] - 3*b_std[:,2], b_mean[:,2] + 3*b_std[:,2], alpha=0.3)
        axs[1, 0].set_ylabel(r'$\eta$')
        axs[1, 1].plot(x_test, b_mean[:,3])
        axs[1,1].fill_between(x_test.ravel(), b_mean[:,3] - 3*b_std[:,3], b_mean[:,3] + 3*b_std[:,3], alpha=0.3)
        axs[1, 1].set_ylabel(r'$Q_{pot} \mathrm{J/kg}$')


    if case == 'homogenization':
        fig, axs = plt.subplots(1, 2,figsize=(8, 4))
        # make the plots tight
        fig.tight_layout(pad=2.0)
        axs[0].plot(x_test, b_mean[:,0])
        axs[0].fill_between(x_test.ravel(), b_mean[:,0] - 2*b_std[:,0], b_mean[:,0] + 2*b_std[:,0], alpha=0.3)
        axs[0].set_ylabel('$E_{paste}$, Pa')
        axs[1].plot(x_test, b_mean[:,1])
        axs[1].fill_between(x_test.ravel(), b_mean[:,1] - 2*b_std[:,1], b_mean[:,1] + 2*b_std[:,1], alpha=0.3)
        axs[1].set_ylabel('$f_{c,paste}$, Pa')

    for ax in axs.flat:
            ax.set(xlabel=r'$r_{sb}$')
            ax.grid()
            # skip if the below if axis is log scale
            if ax.get_yscale() == 'log':
                continue
            else:
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # common legend at the bottom of the plot
    
    if save_path is not None:
        plt.savefig(save_path + 'learnt_prior_predicted_stats_' + case+ datetime+'.pdf')
    plt.show()

    # plot covariance matrix
    #fig, ax = plt.subplots(1,1)
    #ax.imshow(prior_model.cov().detach().numpy(),cmap='hot', interpolation='nearest')
    # colour axis
    # labels for the cov matrixc heatmap
    if case == 'hydration':
        labels = ['$B_1$','$B_2$', r'$\eta$', r'$Q_{pot}$']
    elif case == 'homogenization':
        labels = ['$E_{paste}$','$f_{c,paste}$']
    else:
        NotImplementedError
    sb.heatmap(prior_model.cov().detach().numpy(),annot=True, xticklabels=labels, yticklabels=labels)
    if save_path is not None:
        plt.savefig(save_path + 'covariance_matrix_'+datetime+'.pdf')
    
    plt.show()
   

def prob_hydration_solver_output(NN_model:object,NN_state_dict:str,cov_params:list,latent_dim:int,
                                 temp_key :str = '20C', save_path=None):
    # GET THE PRIOR MODEL

        # load the state dictionary
    NN_model.load_state_dict(th.load(NN_state_dict))

    # covert the values in csv file to a list
    #cov_params = np.genfromtxt(cov_pararms, delimiter=',').tolist()
    prior_model = prior(mean=NN_model, cov_params=cov_params
                           , cov_type='full',latent_dim=latent_dim)
    file_path = 'usecases/optimization_paper/calibration_data/Excel_files/hydration_data_processed.xlsx'
    df = process_hydration_data(file_path)

    # get the optimized values
    x = th.tensor([[0.0],[0.3],[0.5],[0.85]])
    #b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_hydration.npy')

    pred_sample =prior_model.sample(x,n_samples=100)
    

    fig, ax = plt.subplots(1,1)

    ratio_keys = ['CP0','CP30','CP50','CP85']
    inp_solver = {}
    # extrat the first two characters from the string as int
    inp_solver['T_rxn'] = int(temp_key[:2])
    #inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = df[(temp_key,'CP0','Age')]
    hyd_solver = HydrationSolverWrapper()

    Q_mean = []
    Q_std = []
    for i in range(len(ratio_keys)):
        Q_tmp = []
        for j in range(pred_sample.shape[0]):
            Q_tmp.append(hyd_solver.solve(pred_sample[j,i,:],inp_solver))
        # stack to make a 2d array and take the mean along the rows
        Q_mean.append(np.mean(np.vstack(Q_tmp),axis=0))
        Q_std.append(np.std(np.vstack(Q_tmp),axis=0))
        
    colours = ['blue','orange','green','red']
    labels_exp = [r'$\bm{\hat{Q}}_{r_{sb}=0.0}$',r'$\bm{\hat{Q}}_{r_{sb}=0.30}$',r'$\bm{\hat{Q}}_{r_{sb}=0.50}$',r'$\bm{\hat{Q}}_{r_{sb}=0.85}$']
    labels_pred = [r'$\bm{Q}_{r_{sb}=0.0}$',r'$\bm{Q}_{r_{sb}=0.30}$',r'$\bm{Q}_{r_{sb}=0.50}$',r'$\bm{Q}_{r_{sb}=0.85}$']
    for i in range(len(ratio_keys)):
        ax.plot(df[(temp_key,ratio_keys[i],'Age')], df[(temp_key,ratio_keys[i],'Q')],'x', 
                label=labels_exp[i])
        # label with sharp X marker

        ax.plot(df[(temp_key,'CP0','Age')],Q_mean[i],label=labels_pred[i], color=colours[i])
        ax.fill_between(df[(temp_key,'CP0','Age')].ravel(), Q_mean[i] - 2*Q_std[i], Q_mean[i] + 2*Q_std[i], alpha=0.3, color = colours[i])
    ax.legend()
    ax.set_xlabel('Age (s)')
    ax.set_ylabel(r'Cum. Heat of hydration $\bm{Q}$ (J/gh)')
    ax.set_title(r'$T_{rxn}=$'+temp_key[:2]+r'$^{\circ}C$')
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax.grid()
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
        # the legend goes outside the plot, so we need to make the plot a bit wider
    plt.subplots_adjust(right=0.75)
    if save_path is not None:
        plt.savefig(save_path + 'hydration_solver_output_comparison_' + datetime+ '.pdf')
    plt.show()

def prob_homogenization_solver_output(NN_model:object,NN_state_dict:str,cov_params:list,latent_dim:int,save_path=None):
    # load the state dictionary
    NN_model.load_state_dict(th.load(NN_state_dict))

    # covert the values in csv file to a list
    #cov_params = np.genfromtxt(cov_pararms, delimiter=',').tolist()

    prior_model = prior(mean=NN_model, cov_params=cov_params
                           , cov_type='full',latent_dim=latent_dim) # pass the last value of the cov_params
    x_test = th.arange(0,1.0,0.01).reshape(-1,1)
    pred_sample =prior_model.sample(x_test,n_samples=150)

    # load observed data
    path_to_csv = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'
    data_dict = process_homogenization_data(path_to_csv=path_to_csv)
    obs = [np.array(data_dict['E'])*1e06 ,np.array(data_dict['fc(Mpa)'])*1e06 ]
    # scaling back 
    #pred_sample[:,:,0] = pred_sample[:,:,0]*1e09
    #pred_sample[:,:,1] = pred_sample[:,:,1]*1e07

    homogenization_solver = HomogenizationSolverWrapper()
    y_solver = np.zeros((pred_sample.shape[1],pred_sample.shape[2]))
    z_pred_mean = []
    z_pred_std = []
    for i in range(pred_sample.shape[1]):
        tmp = []
        for j in range(pred_sample.shape[0]):
            tmp.append(homogenization_solver.solve(pred_sample[j,i,:]))
        # stack to make a 2d array and take the mean along the rows
        z_pred_mean.append(np.mean(np.vstack(tmp),axis=0))
        z_pred_std.append(np.std(np.vstack(tmp),axis=0))
                       
    # stack the lists
    z_pred_mean = np.vstack(z_pred_mean)
    z_pred_std = np.vstack(z_pred_std)

    x_test = x_test.detach().numpy()
    # plot
    fig, axs = plt.subplots(1, 2,figsize=(8, 4))
    # make the plots tight
    fig.tight_layout(pad=2.5)
    axs[0].plot(x_test, z_pred_mean[:,0])
    axs[0].fill_between(x_test.ravel(), z_pred_mean[:,0] - 2*z_pred_std[:,0], z_pred_mean[:,0] + 2*z_pred_std[:,0], alpha=0.3)
    axs[0].plot(data_dict['x'], obs[0],'x',label='observed')
    axs[0].set_ylabel('$E_{c}$, Pa')
    axs[1].plot(x_test, z_pred_mean[:,1])
    axs[1].fill_between(x_test.ravel(), z_pred_mean[:,1] - 2*z_pred_std[:,1], z_pred_mean[:,1] + 2*z_pred_std[:,1], alpha=0.3)
    axs[1].plot(data_dict['x'], obs[1],'x',label='observed')
    axs[1].set_ylabel('$f_{c}$, Pa')

    for ax in axs.flat:
        ax.grid()
        ax.set(xlabel=r'$r_{sb}$')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    if save_path is not None:
        plt.savefig(save_path + 'homogenization_solver_output_comparison'+datetime+'.pdf')

    plt.show()

if __name__ == '__main__':
    # path_csv = 'cov_parameters2023_08_17-12_03_47_PM.csv'
    # y_label = [r'$\phi_{cov}$']
    # plot_data(path_csv, labels=y_label)

    # path_csv = 'EM_results2023_08_17-12_03_47_PM.csv'
    # labels = [r'-$\mathcal{F}$',r'$||\nabla_{\phi_{nn}}\mathcal{F}||$',r'$||\nabla_{\phi_{cov}}\mathcal{F}||$']
    # plot_data(path_csv, labels=labels)

    nn_model = NN_mean(input_dim=1, output_dim=4, hidden_dim=20)
    #cov_params = 'cov_parameters2023_08_17-12_03_47_PM.csv'
    #nn_state_dict = 'NN_state_dict_till_itr_150.pth'
    #cov_params = np.genfromtxt(cov_params, delimiter=',').tolist()

    nn_state_dict = 'lebedigital/demonstrator_calibration/misc/nn_mean_hydration.pt'
    cov_params = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    viz_learnt_prior_model(nn_model,nn_state_dict,cov_params,latent_dim=4, case='hydration')

    # viz learnt model for homogenization
    nn_model = NN_mean(input_dim=1, output_dim=2, hidden_dim=20)
    nn_state_dict = 'usecases/optimization_paper/model_learning/homogenization/exp_2/NN_state_dict_till_itr_150_2023_08_22-07_15_51_PM.pth'
    cov_path = 'usecases/optimization_paper/model_learning/homogenization/exp_2/cov_parameters2023_08_22-07_15_51_PM.csv'
    cov_params = np.genfromtxt(cov_path, delimiter=',').tolist()[-1]
    viz_learnt_prior_model(nn_model,nn_state_dict,cov_params,latent_dim=2,case='homogenization',save_path='lebedigital/demonstrator_calibration/misc/')
    
    
    prob_homogenization_solver_output(nn_model,nn_state_dict,cov_params,latent_dim=2,save_path='lebedigital/demonstrator_calibration/misc/')
    #prob_hydration_solver_output(nn_model,nn_state_dict,cov_params,latent_dim=4,save_path='learnt_prior_')





