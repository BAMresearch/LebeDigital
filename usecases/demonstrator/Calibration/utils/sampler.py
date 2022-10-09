import numpy as np
from tqdm import tqdm
import scipy.stats as ss
# TODO: implement component wise https://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
#https://utstat.toronto.edu/craiu/Talks/uqam_talk.pdf
class random_walk_metropolis:
    def __init__(self,target_logprob):
        self._target_log_prob = target_logprob
        self.scale_cov = None
        self.acceptance_ratio = None


    def run(self,N, cov_proposal, x0,burnin =None, **kwargs):
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
        assert cov_proposal.ndim == 2, "Full cov matrix must be supplied"
        dimx = np.size(x0)
        logp = self._target_log_prob(x0,**kwargs)
        accepted = 0

        X_chain = np.zeros((N, dimx))
        scale_cov = 1. # start with no proposal cov scaling
        for n in tqdm(range(N)):
            # The proposal distribution goes here
            # x_proposed = x + stepsize*np.random.normal(0,1,dimx)
            cov_proposal_scaled = scale_cov*cov_proposal
            x_proposed = ss.multivariate_normal(mean=x,cov=cov_proposal_scaled).rvs()
            logp_proposed = self._target_log_prob(x_proposed,**kwargs)   # Target density

            #if np.random.uniform() <= logp_proposed/logp:  #(as we took log of the acceptance ratio)
            #alpha = min(1, np.exp(logp_proposed - logp))
            #if np.random.rand() <= alpha: # accept
            if np.log(np.random.uniform())<= logp_proposed - logp:
                # accept
                x=x_proposed
                logp = logp_proposed
                accepted += 1
            #if (n>1):#and (n%20 == 0)): # to avoid division by 0
            #    scale_cov = self._tune_scale_covariance(scale_covariance=scale_cov,accept_rate=accepted/n)
            X_chain[n,:] = x
            self.scale_cov = scale_cov
        self.acceptance_ratio = accepted/N
        print("Acceptance ratio: {} and cov scale: {}".format(self.acceptance_ratio, scale_cov))
        if burnin is not None:
            #TODO acceptance rate shouldnt include burnin samples
            X_chain = X_chain[burnin:,:]
        return X_chain

    def _tune_scale_covariance(self,scale_covariance, accept_rate):
        """
        Tune the acceptance rate according to the last tuning interval. If higher acceptance rate , means
        you need to expand you search field or increase variance(its too small currently)

        The goal is an acceptance rate within 20\% - 50\%.
        The (acceptance) rate is adapted according to the following rule:

            Acceptance Rate    Variance adaptation factor
            ---------------    --------------------------
            <0.001                       x 0.1
            <0.05                        x 0.5
            <0.2                         x 0.9
            >0.5                         x 1.1
            >0.75                        x 2
            >0.95                        x 10

        The implementation is modified from [1].

        Reference:
        [1]: https://github.com/pymc-devs/pymc3/blob/master/pymc3/step_methods/metropolis.py
        """
        if accept_rate < 0.001:
            scale_covariance = 0.1*scale_covariance
        if ((accept_rate>=0.001) and (accept_rate<0.05)):
            scale_covariance = 0.5*scale_covariance
        if ((accept_rate >= 0.05) and (accept_rate < 0.2)):
            scale_covariance = 0.9 * scale_covariance
        if ((accept_rate >= 0.5) and (accept_rate < 0.75)):
            scale_covariance = 1.1 * scale_covariance
        if ((accept_rate >= 0.75) and (accept_rate < 0.95)):
            scale_covariance = 2 * scale_covariance
        if (accept_rate >= 0.95):
            scale_covariance = 10 * scale_covariance



        # scale_covariance = np.where(accept_rate < 0.001, scale_covariance * 0.1, scale_covariance)
        # scale_covariance = np.where(
        #     (accept_rate >= 0.001) * (accept_rate < 0.05), scale_covariance * 0.5, scale_covariance
        # )
        # scale_covariance = np.where(
        #     (accept_rate >= 0.05) * (accept_rate < 0.2), scale_covariance * 0.9, scale_covariance
        # )
        # scale_covariance = np.where(
        #     (accept_rate > 0.5) * (accept_rate <= 0.75), scale_covariance * 1.1, scale_covariance
        # )
        # scale_covariance = np.where(
        #     (accept_rate > 0.75) * (accept_rate <= 0.95), scale_covariance * 2.0, scale_covariance
        # )
        # scale_covariance = np.where((accept_rate > 0.95), scale_covariance * 10.0, scale_covariance)

        return scale_covariance