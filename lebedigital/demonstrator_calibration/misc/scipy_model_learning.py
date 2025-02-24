# %%
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper,HomogenizationSolverWrapper
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data, process_homogenization_data
from lebedigital.demonstrator_calibration.parametric_model import train_NN, NN_mean
import numpy as np
import torch as th
from scipy.optimize import minimize
from scipy import odr
import os
from matplotlib import pyplot as plt
# use latex with matplotlib
plt.rc('text', usetex=True)
from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")


# set deafult data type to float32
th.set_default_dtype(th.float64)

#%%
def loss(params, inp_obs:dict,forward_model:callable, obs:np.array):

    z_pred = forward_model(latents=params, inp_solver=inp_obs)
    # normalisation
    #Q_pred = (Q_pred- np.mean(Q_pred))/(np.std(Q_pred) + 1e-07)
    #Q_exp = (Q_exp- np.mean(Q_exp))/(np.std(Q_exp) + 1e-07)
    #assert Q_exp.shape == Q_pred.shape
    # convert obs and z_pred to array if they are not
    obs = np.array(obs)
    z_pred = np.array(z_pred)
    obj = np.sqrt(np.mean((z_pred - obs) ** 2)) # RMS
    # plot z_pred and obs to see if they are aligned
    # plt.plot(z_pred, label='z_pred')
    # plt.plot(obs, label='obs')
    # plt.legend()
    # plt.show()
    return obj

def optimize_for_all(b_init,solver:callable,data:dict, ratio_keys):
    b_opt = []
    for i,v in enumerate(ratio_keys):
        inp_solver = {}
        inp_solver['T_rxn'] = 20
        inp_solver['time_list'] = data[('20C',v,'Age')] # temp. hardcoded to the hydration model
        obs = data[('20C',v,'Q')]
        res = minimize(loss, x0 = b_init ,args=(inp_solver,solver,obs), method='Nelder-Mead', options = {'maxiter': 500} )
        #x_init = res.x
        b_opt.append(res.x)
    b_opt = np.stack(b_opt)
    return b_opt

def optimize_for_all_homogenization(b_init, solver:callable, data:dict):
    b_opt = []
    for i,v in enumerate(data['x']):
        inp_solver = None
        obs = [data['E'][i]*1e06 ,data['fc(Mpa)'][i]*1e06 ]# convert to Pa from Mpa
        res = minimize(loss, x0 = b_init ,args=(inp_solver,solver,obs), method='Nelder-Mead', options = {'maxiter': 500} )
        #x_init = res.x
        b_opt.append(res.x)
    b_opt = np.stack(b_opt)
    return b_opt

if __name__ == '__main__':
    #%%
    #file_location = os.path.dirname(os.path.realpath(__file__))
    #path = file_location + '/Excel_files/hydration_data_processed.xlsx'
    optimize = False
    train = True
    plot_solver_output = False
    optimize_homogenization = False
    plot_homogenization_output = False
    scipy_poly_regression_homogenization = False
    nn_fit_homogenization = False


    if optimize_homogenization:
        path_to_csv = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'

        data_dict = process_homogenization_data(path_to_csv=path_to_csv)
        homogenization_solver = HomogenizationSolverWrapper()

        b_init = [30e9,30e6]
        b_opt = optimize_for_all_homogenization(b_init, homogenization_solver.solve, data_dict)
        np.save('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_'+datetime+'.npy',b_opt)

        fig, axs = plt.subplots(1, 2)
        fig.tight_layout(pad=3.0)
        axs[0].plot(data_dict['x'], b_opt[:,0], '-*')
        axs[0].set_title('E')
        axs[1].plot(data_dict['x'],b_opt[:,1], '-*')
        axs[1].set_title('fc')
        plt.savefig('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_'+datetime+'.png')
        plt.show()

    if plot_homogenization_output:
        path_to_csv = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'
        data_dict = process_homogenization_data(path_to_csv=path_to_csv)
        homogenization_solver = HomogenizationSolverWrapper()
        inp_solver = None

        b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_2023_08_21-05_47_52_PM.npy')
        
        fig, ax = plt.subplots(1, 2)
        fig.tight_layout(pad=3.0)
        # make the subplots square
        ax[0].set_aspect('equal', 'box')
        ax[1].set_aspect('equal', 'box')
        for i,v in enumerate(data_dict['x']):
            Z_pred = homogenization_solver.solve(latents=b_opt[i,:], inp_solver=None)
            obs = [data_dict['E'][i]*1e06 ,data_dict['fc(Mpa)'][i]*1e06 ]# convert to Pa from Mpa

            ax[0].plot(Z_pred[0],obs[0],'x-',color='blue')
            ax[0].set_title('E')
            ax[0].set_ylabel('observed')
            ax[0].set_xlabel('predicted')
            ax[1].plot(Z_pred[1],obs[1],'x-',color='blue')
            ax[1].set_title('fc')            
            ax[1].set_ylabel('observed')
            ax[1].set_xlabel('predicted')
        plt.savefig('lebedigital/demonstrator_calibration/misc/homogenization_prediction_comparission'+datetime+'.png')
        plt.show()

    if scipy_poly_regression_homogenization:
        path_to_csv = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'
        b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_2023_08_21-05_47_52_PM.npy')
        data_dict = process_homogenization_data(path_to_csv=path_to_csv)

        poly_model = odr.polynomial(3)
        data = odr.Data(data_dict['x'],b_opt[:,0])
        odr_obj = odr.ODR(data,poly_model)
        output = odr_obj.run()
        poly = np.poly1d(output.beta[::-1])
        poly_y = poly(data_dict['x'])
        plt.plot(data_dict['x'], poly_y, label="polynomial ODR")

    if nn_fit_homogenization:
        path_to_csv = 'usecases/optimization_paper/calibration_data/Excel_files/homogenization_data_processed_E.csv'
        b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_homogenization_2023_08_21-05_47_52_PM.npy')
        data_dict = process_homogenization_data(path_to_csv=path_to_csv)

        x_train = th.tensor(data_dict['x'])
        # make the above 2d tensor
        x_train = x_train.reshape(-1,1)
        y_train = th.tensor(b_opt)
        # scale the y_train values for training by diving with the magnitude
        y_train[:,0] = y_train[:,0]/1e09
        y_train[:,1] = y_train[:,1]/1e07


        nn_mean = train_NN(NN_mean,x=x_train, y=y_train, epochs=4000, lr=1e-2, hidden_dim=30)

        # generate 0 to 1 qith 0.1 step and do prediction for that
        x_test = th.arange(0,1.0,0.01).reshape(-1,1)
        y_test = nn_mean(x_test)

        # scale back the values
        y_test[:,0] = y_test[:,0]*1e09
        y_test[:,1] = y_test[:,1]*1e07
        y_train[:,0] = y_train[:,0]*1e09
        y_train[:,1] = y_train[:,1]*1e07

        # plot the results column wise
        # covert to numpy
        x_test = x_test.detach().numpy()
        y_test = y_test.detach().numpy()
        x_train = x_train.detach().numpy()
        y_train = y_train.detach().numpy()

        # plot
        fig, axs = plt.subplots(2, 1)
        # make the plots tight
        fig.tight_layout(pad=3.0)
        #axs[0].set_aspect('equal', 'box')
        #axs[1].set_aspect('equal', 'box')

        axs[0].plot(x_test.ravel(), y_test[:,0])
        axs[0].plot(x_train.ravel(), y_train[:,0], 'x')
        axs[0].set_title('E_paste')
        axs[1].plot(x_test.ravel(), y_test[:,1])
        axs[1].plot(x_train.ravel(), y_train[:,1], 'x')
        axs[1].set_title('fc_paste')
        plt.savefig('lebedigital/demonstrator_calibration/misc/nn_mean_homogenization_prediction'+datetime+'.png')
        plt.show()

        homogenization_solver = HomogenizationSolverWrapper()
        # create a new variable y_solver same size as y_test
        y_solver = np.zeros(y_test.shape)
        for i in range(y_test.shape[0]):
            y_solver[i,:] = homogenization_solver.solve(latents=y_test[i,:], inp_solver=None)
        obs = [np.array(data_dict['E'])*1e06 ,np.array(data_dict['fc(Mpa)'])*1e06 ]# convert to Pa from Mpa

        # plot 
        fig, axs = plt.subplots(2, 1)
        # make the plots tight
        fig.tight_layout(pad=3.0)

        axs[0].plot(x_test.ravel(), y_solver[:,0])
        axs[0].plot(data_dict['x'], obs[0], 'x')
        axs[0].set_title('E_concrete')
        axs[1].plot(x_test.ravel(), y_solver[:,1])
        axs[1].plot(data_dict['x'], obs[1], 'x')
        axs[1].set_title('fc_concrete')
        plt.savefig('lebedigital/demonstrator_calibration/misc/nn_mean_homogenization_prediction_concrete'+datetime+'.png')
        plt.show()




    if optimize:
        file_path = 'usecases/optimization_paper/calibration_data/Excel_files/hydration_data_processed.xlsx'
        hydration_data = process_hydration_data(file_path)
        print(hydration_data.keys())

        hydration_solver = HydrationSolverWrapper()
        # inp_solver = {}
        # inp_solver['T_rxn'] = 20
        # inp_solver['time_list'] = hydration_data[('20C','CP0','Age')]
        # q_obs = hydration_data[('20C','CP0','Q')]
        # q_pred = hydration_solver.solve([2.916,2.4229,5.554,5],inp_solver)

        # obj = loss(params=[2.916,2.4229,5.554,5],inp_obs=inp_solver,forward_model=hydration_solver.solve,obs=q_obs)
        
        #b_init = np.array([2.916,2.4229,5.554,5])
        b_init = np.array([2.916E-4,0.0024229,5.554,500e3])

        # log-traansform the parameters
        b_init = np.log(b_init)

        ratio_keys = ['CP0','CP30','CP50','CP85']
        b_opt = optimize_for_all(b_init,hydration_solver.solve,hydration_data,ratio_keys)
        # save b_opt as np.save here in the currente directory
        np.save('lebedigital/demonstrator_calibration/misc/optimization_results_hydration_'+datetime+'.npy',b_opt)

        # subplot for the results
        fig, axs = plt.subplots(2, 2)
        fig.tight_layout(pad=3.0)
        # axs[0, 0].semilogy(b_opt[:,0], '-*')
        # axs[0, 0].set_title('$B_1$')
        # axs[0, 1].semilogy(b_opt[:,1], '-*')
        # axs[0, 1].set_title('$B_2$')
        # axs[1, 0].semilogy(b_opt[:,2], '-*')
        # axs[1, 0].set_title(r'$\eta$')
        # axs[1, 1].semilogy(b_opt[:,3], '-*')
        # axs[1, 1].set_title(r'$Q_{pot}$')
        axs[0, 0].plot(b_opt[:,0], '-*')
        axs[0, 0].set_title('$B_1$')
        axs[0, 1].plot(b_opt[:,1], '-*')
        axs[0, 1].set_title('$B_2$')
        axs[1, 0].plot(b_opt[:,2], '-*')
        axs[1, 0].set_title(r'$\eta$')
        axs[1, 1].plot(b_opt[:,3], '-*')
        axs[1, 1].set_title(r'$Q_{pot}$')
        for ax in axs.flat:
            ax.set(xlabel='slag/cement ratio', ylabel='value')
        plt.savefig('lebedigital/demonstrator_calibration/misc/optimization_results_scipy_hydration_'+datetime+'.png')
        plt.show()
    if plot_solver_output:
        file_path = 'usecases/optimization_paper/calibration_data/Excel_files/hydration_data_processed.xlsx'
        df = process_hydration_data(file_path)

        # get the optimized values
        x = th.tensor([[0.0],[0.3],[0.5],[0.85]])
        #b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_hydration.npy')
        b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_hydration_2023_08_18-12_58_47_PM.npy')
        fig, ax = plt.subplots(1,1)

        ratio_keys = ['CP0','CP30','CP50','CP85']
        inp_solver = {}
        inp_solver['T_rxn'] = 20
        inp_solver['time_list'] = df[('20C','CP0','Age')]
        hyd_solver = HydrationSolverWrapper()


        for i in range(len(ratio_keys)):
            ax.plot(df[('20C',ratio_keys[i],'Age')], df[('20C',ratio_keys[i],'Q')],'+-', 
                    label=ratio_keys[i]+'exp')
            # label with sharp X marker

            ax.plot(df[('20C','CP0','Age')],hyd_solver.solve(b_opt[i,:],inp_solver),'x-',label=ratio_keys[i]+'pred')
        ax.legend()
        ax.set_xlabel('Age (s)')
        ax.set_ylabel('Cum. Heat of hydration (J/gh)')
        fig.savefig('lebedigital/demonstrator_calibration/misc/hydration_solver_output_comparison_'+datetime+'.png')
        plt.show()


    if train:
        x = th.tensor([[0.0],[0.3],[0.5],[0.85]])
        b_tmp = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_hydration.npy')
        y = th.tensor(b_tmp)

        # do a test/tain split. In train, remove the 3rd row and keep all others
        # in test, keep only the 3rd row
        x_train = x[[0,1,2,3],:]
        y_train = y[[0,1,2,3],:]

        # standardize the training data
        
        #y_train = (y_train - th.mean(y_train,dim=0))/th.std(y_train,dim=0)
        y_train[:,1] = y_train[:,1]*1e-03
        y_train[:,1] = th.log(y_train[:,1]) # log transform the B2

        nn_mean = train_NN(NN_mean,x=x_train, y=y_train, epochs=2000, lr=1e-02, hidden_dim=20)

        # save the model
        #th.save(nn_mean.state_dict(), 'lebedigital/demonstrator_calibration/misc/nn_mean_hydration.pt')
        
        # load the model
        #nn_mean = NN_mean(input_dim=x_train.shape[1], output_dim=y_train.shape[1], hidden_dim=30)
        #nn_mean.load_state_dict(th.load('lebedigital/demonstrator_calibration/misc/nn_mean_hydration.pt'))

        # test the model
        y_pred = nn_mean(x[[3],:])
        print(y_pred)

        # generate 0 to 1 qith 0.1 step and do prediction for that
        x_test = th.arange(0,1.0,0.01).reshape(-1,1)
        y_test = nn_mean(x_test)

        # plot the results column wise
        # covert to numpy
        x_test = x_test.detach().numpy()
        y_test = y_test.detach().numpy()
        x_train = x_train.detach().numpy()
        y_train = y_train.detach().numpy()

        fig, axs = plt.subplots(2, 2)
        # make the plots tight
        fig.tight_layout(pad=3.0)
        axs[0, 0].plot(x_test, y_test[:,0])
        axs[0,0].plot(x_train, y_train[:,0], 'x')
        axs[0, 0].set_title('$B_1$')
        axs[0, 1].plot(x_test, y_test[:,1])
        axs[0,1].plot(x_train, y_train[:,1], 'x')
        #axs[0, 1].set_title('$B_2$')
        axs[0, 1].set_title('$\log(B_2)$')
        axs[1, 0].plot(x_test, y_test[:,2])
        axs[1,0].plot(x_train, y_train[:,2], 'x')
        axs[1, 0].set_title(r'$\eta$')
        axs[1, 1].plot(x_test, y_test[:,3])
        axs[1,1].plot(x_train, y_train[:,3], 'x')
        axs[1, 1].set_title(r'$Q_{pot}$')
        for ax in axs.flat:
            ax.set(xlabel='slag/cement ratio', ylabel='value')
        plt.savefig('lebedigital/demonstrator_calibration/misc/nn_mean_hydration_prediction.png')
        plt.show()




