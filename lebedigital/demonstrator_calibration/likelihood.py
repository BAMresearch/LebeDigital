import torch as th
import numpy as np
from abc import ABC, abstractmethod
from lebedigital.demonstrator_calibration.forward_solvers import HydrationSolverWrapper

class likelihood(ABC):
    """Base class for likelihoods
    """
    def __init__(self,forward_model:callable,sigma:float):
        """_summary_

        Parameters
        ----------
        forward_model : callable
            physics based model which acts as an observation operator. accepts kwargs as inputs.
            needs 'latents' and 'inp_solver' as inputs. check forwardS_solvers.py for more info
        sigma : _type_
            s.d of the likelihood. hardcoded to be a scalar for now
        """
        self.forward_model = forward_model
        self.sigma = sigma

    @abstractmethod
    def log_prob(self,observed,**kwargs):
        pass

class gaussian_likelihood:
    def __init__(self,forward_model:callable,sigma):
        """_summary_

        Parameters
        ----------
        forward_model : callable
            physics based model which acts as an observation operator. accepts kwargs as inputs
        sigma : _type_
            s.d of the likelihood. hardcoded to be a scalar for now
        """
        self.forward_model = forward_model
        self.sigma = sigma
    
    def log_prob(self,observed,**kwargs):
        z_pred = self.forward_model(**kwargs) # forward model needs 'latents' and 'inp_solver' as inputs
        # TODO: infer sigma also
        #cov = th.eye(z_pred.shape[0])*self.sigma**2
        #if sigma is a list, create a diagonal covariance
        if isinstance(self.sigma,list):
           cov = th.diag(th.tensor(self.sigma)**2)
        else:
            cov = th.eye(z_pred.shape[0])*self.sigma**2
        dist = th.distributions.MultivariateNormal(th.tensor(z_pred),cov)
        return dist.log_prob(th.tensor(observed))

if __name__ == '__main__':
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = [0,5000,10000,20000,100000]

    # -- latents -----
    b = [2.916,2.4229,5.554,5]
    hydration_solver = HydrationSolverWrapper()
    obs =[  0.        ,   2.14549938,   7.1823244 ,  34.34254352,
       233.33527714]
    # multiple every element of obs by 2
    obs = [i*2 for i in obs]
    
    lkl = gaussian_likelihood(forward_model=hydration_solver.solve,sigma=1)
    lp = lkl.log_prob(observed=obs,latents=b,inp_solver=inp_solver)
    print(f'the log_prob is {lp}')