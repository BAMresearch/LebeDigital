
import torch as th
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import copy
import pandas as pd
import csv

from datetime import datetime
datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
import matplotlib as mpl
from matplotlib import rc

from lebedigital.demonstrator_calibration.parametric_model import NN_mean, train_NN
from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.demonstrator_calibration.likelihood import gaussian_likelihood
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper
from usecases.optimization_paper.calibration_data.data_handling import process_hydration_data
from lebedigital.demonstrator_calibration.sampler import MCMC_DRAM

# set torch deafult data type to float32

class VBEM:
    """class implementing the Variational Bayes Expectation Maximization algorithm. 
    Ugly copy of the VBEM class in VBEM.py for the homogenization."""
    def __init__(self,prior:prior,forward_model:callable,likelihood:gaussian_likelihood, 
                 model_prior_mean:NN_mean, prior_cov_params:list, sigma_likelihood:float, latent_dim:int, 
                 dataframe_observed_data:dict,no_observed_data_pair:int,b_init:list,pre_train:bool=True,lr=1e-2):
        
        #assert len(b_init) == latent_dim, 'the length of the latents must be equal to the latent dim'
        self.b_init = b_init
        self.df = dataframe_observed_data
        # initialize the Neural Nets
        self.NN_mean = model_prior_mean

        # initialize the parameters of the NN model with pretraining
        if pre_train:
            self.nn_mean = self._pre_train_network()
        else:
            self.nn_mean = model_prior_mean
        
        # initialize the prior and likelihood
        self.latent_dim = latent_dim
        self.prior_cov_params = th.tensor(prior_cov_params,requires_grad=True)
        self.prior = prior(mean=self.nn_mean, cov_params=self.prior_cov_params
                           , cov_type='full',latent_dim=self.latent_dim)

        self.forward_model = forward_model
        self.sigma_likelihood = sigma_likelihood
        self.likelihood = likelihood(self.forward_model, self.sigma_likelihood)

        # define the optimizer
        # add learning rate schedular

        self.optimizer = th.optim.Adam([{'params': self.prior.mean.parameters()},
                                        {'params': self.prior.para_cov_torch}], lr=lr,weight_decay=1e-03)
        #self.scheduler = th.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=1)
        
        # # assert if the observed_data has certain keys 'x and z'
        # assert 'x' in observed_data.keys(), 'the observed data must have the keys x and z'
        # assert 'z' in observed_data.keys(), 'the observed data must have the keys x and z'
  
        # # assert the values of the keys are are atleast 2d arrays
        # assert len(observed_data['x'].shape) == 2, 'the observed data input must be atleast 2d arrays'
        # assert len(observed_data['z'].shape) == 2, 'the observed data input must be atleast 2d arrays'
        # self.observed_data = observed_data
        
        self.no_obs = no_observed_data_pair
        # tmp variables
        self.x_tmp = None
        self.z_tmp = None
        self.inp_solver_tmp = None
        self.x = np.array(dataframe_observed_data['x']).reshape(-1,1).tolist() # the input data
        #self.x = [[0.3]]

        # define the data holders
        self.q_b_list = []
        self.grad_prior_parameters_holder = []
        self.loss_holder = []
    
    def _pre_train_network(self):
        #TODO this is hugly hard coded, fix it
        "pre train the network for better weight initialization"
            # run the pretraining for 1 dim input and 4 dim output synthetic data 
        #x = th.tensor([[0.3],[0.6]])
        x = th.tensor(self.df['x']).reshape(-1,1)
        #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
        #y = th.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
        y = th.tensor(self.b_init)
        y = y + 0.05*th.randn_like(y)
        nn_mean = train_NN(self.NN_mean,x, y, epochs=700, lr=1e-2, hidden_dim=20)
        return nn_mean
    
    def _posterior_model(self,b:list):
        """ b is the list of latents"""
        assert len(b) == self.latent_dim, 'the length of the latents must be equal to the latent dim'
        assert self.z_tmp is not None, 'the observed data must be set before calling this function'
        #assert self.inp_solver_tmp is not None, 'the input solver must be set before calling this function'
        assert self.x_tmp is not None, 'the input data must be set before calling this function'

        return self.prior.log_prob(x=self.x_tmp,b=b) + self.likelihood.log_prob(observed=self.z_tmp,latents = b,
                                                               inp_solver=self.inp_solver_tmp)
    def _temp_input_hydration_model(self,x:int):
        """set the temporary input for the homogenization model"""
        #TODO  do the unit conversion outside of this function
        self.z_tmp = [self.df['E'][x]*1e06, self.df['fc(Mpa)'][x]*1e06] # converting from Mpa to Pa
        self.x_tmp = self.x[x]



    def _E_step(self, no_samples):
        """run the E step of the VBEM algorithm"""

        for i in range(self.no_obs):
            # parallelize the loop using openMPI:
            # https://mpi4py.readthedocs.io/en/stable/tutorial.html

            # initlialize the values for log_posterior
            self._temp_input_hydration_model(i)
            samples_df = MCMC_DRAM(log_func=self._posterior_model,n_dim=self.latent_dim, no_samples=no_samples,
                                   x_init=self.b_init[i])
            #samples_df = MCMC_DRAM(log_func=self._posterior_model,n_dim=self.latent_dim, no_samples=100)

            # covert to 2D array
            q_b = samples_df.to_numpy()[:,1:]
            
            # append to the list
            self.q_b_list.append(q_b)

            # update the b_init
            self.b_init[i] = q_b[-1,:].tolist()
        q_b_list = self.q_b_list
        # clear the list
        self.q_b_list = []
        return q_b_list

    def M_step(self, no_steps:int, no_samples:int=100,q_sample_test=None):
        # TODO: write a test with identical samples of a specific latent from scipy-opt
        # and check if the NN is able to recover the same latent.
        for i in range(no_steps):
            self.optimizer.zero_grad()
            # run the E step and collect a list of latents
            if q_sample_test is not None:
                #q_b_samples = [np.array([[1.7, 1.43, 1.56, 1.8],[1.7, 1.43, 1.56, 1.8],[1.7, 1.43, 1.56, 1.8]])]
                q_b_samples = q_sample_test
            else:
                if i==0:
                    q_b_samples = self._E_step(no_samples=5*no_samples) # for good initialization
                else:
                    q_b_samples = self._E_step(no_samples=no_samples)
            print(f'q_b_samples:{q_b_samples}, q length:{len(q_b_samples)}, q_shape: {q_b_samples[0].shape}')
            # get the grad estimate for the latent parameters
            #grad_NN, grad_cov = self.prior.grad_estimate_score_function(self.observed_data['x'],q_b_samples,return_grad=True)
            E_log_prob = self.prior.grad_estimate_score_function(self.x,q_b_samples,return_grad=False)
            
            # run the optimizer step
            obj = -E_log_prob
            obj.backward()
            self.optimizer.step()
            #self.scheduler.step()

            with th.no_grad():
                grad_norm_nn = th.norm(th.cat([p.grad.flatten() for p in self.prior.mean.parameters()]))
                grad_norm_cov = th.norm(self.prior.para_cov_torch.grad.flatten())
                nn_parameters = th.cat([p.flatten() for p in self.prior.mean.parameters()]).detach().numpy()
                cov_parameters = self.prior.para_cov_torch.detach().numpy()
            if i % 1 == 0:
                print(f'iteration {i}, objective : {obj}, grad norm: {grad_norm_nn}, cov: {self.prior.para_cov_torch}, \
                      cov_grad: {self.prior.para_cov_torch.grad}')
                #print(f'cov params: {self.prior.para_cov_torch.grad}')

            # open a .csv file and write the results for each iteration
            with th.no_grad():
                #path = 
                with open('EM_results'+datetime+ '.csv', 'a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow([obj.item(),grad_norm_nn.item(),grad_norm_cov.item()])
                #pd.DataFrame(nn_parameters).to_csv('NN_parameters'+datetime+'.csv',mode='a',header=False)
                #pd.DataFrame(cov_parameters).to_csv('cov_parameters'+datetime+'.csv',mode='a',header=False)
                # append the NN wreights to dataframe and save to csv file row wise
                df = pd.DataFrame(nn_parameters).T
                df.to_csv('NN_parameters'+datetime+'.csv',mode='a',header=False,index=False)
                # append the cov parameters to dataframe and save to csv file row wise
                df = pd.DataFrame(cov_parameters).T
                df.to_csv('cov_parameters'+datetime+'.csv',mode='a',header=False,index=False)

            # saving state_dict of the model and optimizer, can be used to resume training 
            if i % 50 == 0:
                th.save(self.prior.mean.state_dict(), 'NN_state_dict_till_itr_'+ str(i) +'_'+datetime+'.pth')
                th.save(self.optimizer.state_dict(), 'optimizer_state_dict_till_itr_'+ str(i) +'_'+datetime+'.pth')
        
        # script the NN to call without instantiating the class
        model_scripted = th.jit.script(self.prior.mean)
        model_scripted.save('NN_mean_scripted_'+datetime+'.pt')



#%%
# ----------------------------

if __name__ == '__main__':

    hydration_solver = HydrationSolverWrapper()
    cov_param = th.tensor([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],requires_grad=True)
    #cov_param = th.tensor([10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0],requires_grad=True)
    observation_data = {'x':np.array([[0.3]]),
                        'z':np.random.normal(size=(2,4))} # z will take the solver output
    data_path = 'usecases/optimization_paper/calibration_data/Excel_files/hydration_data_processed.xlsx'
    df_data = process_hydration_data(data_path)
    b_opt = np.load('lebedigital/demonstrator_calibration/misc/optimization_results_hydration.npy')
    b_init = b_opt.tolist()
    #b_init[3] = b_opt.tolist()[2]
    # TODO: 1. b's need to be normalized.
    # TODO: the forth data point looks fishy.

    vbem = VBEM(prior=prior,forward_model=hydration_solver.solve,likelihood=gaussian_likelihood,
                model_prior_mean=NN_mean,prior_cov_params=cov_param, sigma_likelihood=1, latent_dim=4,
                dataframe_observed_data=df_data,no_observed_data_pair=4,b_init=b_init,
                pre_train=True)
    #vbem.M_step(5)
    
    # save the trained NN in a file
    #th.save(vbem.prior.mean.state_dict(), 'NN_mean.pt')
    q_sample_test = [np.array([[1.7, 1.43, 1.56, 1.8],[1.7, 1.43, 1.56, 1.8],[1.7, 1.43, 1.56, 1.8]])]
    vbem.M_step(1000, q_sample_test=q_sample_test)

    pred = vbem.prior.mean(th.tensor([[0.3]]))
    print(pred)
    pred_sample = vbem.prior.sample(th.tensor([[0.3]]),n_samples=100)
    # get the mean and variance of the samples
    print(f'The sample mean is {np.mean(pred_sample,axis=0)}') # [[1.70475122 1.42880646 1.57344777 1.79734188]]
    print(f'the sample var is {np.var(pred_sample,axis=0)}') # [[0.0019244  0.00191447 0.04464581 0.01586459]]

    print('done')


    # ----------------------------
    # run the pretraining for 1 dim input and 4 dim output synthetic data 
    # x = th.tensor([[0.3],[0.6]])
    # #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
    # y = th.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
    # nn_mean = train_NN(NN_mean,x, y, epochs=2000, lr=1e-2, hidden_dim=10)

    # temp_cov_param = th.tensor([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],requires_grad=True)
    # optimizer = th.optim.Adam([{'params': nn_mean.parameters()},{'params': temp_cov_param}], lr=1e-2)

    # prior_ = prior(nn_mean, cov_params=temp_cov_param, cov_type='full',latent_dim=4)

    # # values at which the gradient is expected to be almost close to zero
    # for i in range(5):
    #     b_list_new = [np.array([[2.7, 2.43, 5.56, 4.8],[2.7, 2.43, 5.56, 4.8],[2.7, 2.43, 5.56, 4.8]])]
    #     grad_direct_mean, grad_direct_cov = prior_.grad_estimate_score_function([[0.6]],b_list_new,return_grad=True)
    #     optimizer.zero_grad()
    #     # pass nn_mean.parameters() and temp_cov_param to the optimizer

    #     print('parameters bedfore the step')
    #     for i,p in enumerate(nn_mean.parameters()):
    #         print(p)
    #     print(f'the cov parameters : {temp_cov_param}')
    #     for i,p in enumerate(nn_mean.parameters()):
    #         p.grad = grad_direct_mean[i]
    #     #temp_cov_param.grad = grad_direct_cov
    #     optimizer.step()
    #     print('parameters after the step')
    #     for i,p in enumerate(nn_mean.parameters()):
    #         print(p)
    #     print(f'the cov parameters : {temp_cov_param}')

    #     print('the gradient of the parameters')
    #     for p in nn_mean.parameters():
    #         print(p.grad, temp_cov_param.grad)