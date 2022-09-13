import numpy as np
from tqdm import tqdm
import scipy.stats as ss

class random_walk_metropolis:
    def __init__(self,target_logprob):
        self._target_log_prob = target_logprob

    def run(self,N, stepsize, x0,burnin =None, **kwargs):
        """

        Parameters
        ----------
        N :
        stepsize :
        x0 :
        obs_data :
        i : Index of the observed datapair
        burnin :

        Returns
        -------

        """
        x = x0  #Intial value for mut essentially/start with {0}
        dimx = np.size(x0)
        logp = self._target_log_prob(x0,**kwargs)
        accepted = 0

        X_chain = np.zeros((N, dimx))

        for n in tqdm(range(N)):
            # The proposal distribution goes here
            # x_proposed = x + stepsize*np.random.normal(0,1,dimx)
            x_proposed = ss.multivariate_normal(mean=x,cov=np.diag(stepsize)).rvs()
            logp_proposed = self._target_log_prob(x_proposed,**kwargs)   # Target density

            #if np.random.uniform() <= logp_proposed/logp:  #(as we took log of the acceptance ratio)
            #alpha = min(1, np.exp(logp_proposed - logp))
            #if np.random.rand() <= alpha: # accept
            if np.log(np.random.uniform())<= logp_proposed - logp:
                # accept
                x=x_proposed
                logp = logp_proposed
                accepted += 1
            X_chain[n,:] = x

        print("Acceptance ratio: {}".format(accepted / N))
        if burnin is not None:
            X_chain = X_chain[burnin:,:]
        return X_chain