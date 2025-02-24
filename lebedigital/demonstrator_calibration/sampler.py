#%%
import numpy as np
from tqdm import tqdm
import scipy.stats as ss
import paramonte as pm
import matplotlib.pyplot as plt
import seaborn as sns

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

# ----------------
#%%
def MCMC_DRAM(log_func:callable, n_dim:int,no_samples=1000,x_init:list =None, **kwargs):
    """
    MCMC using DRAM (Delayed Rejection Adaptive Metropolis)
    https://www.cdslab.org/paramonte/notes/overview/preface/#what-is-paramonte
    Parameters
    ----------
    log_func : callable
        The first agrumnet should be the RVs and the rest are the parameters which should be provided.
    n_dim : int
        _description_
    no_samples : int
        _description_
    x_init : list
        The initial point of the Markov Chain
    """

    pmpd = pm.ParaDRAM()
    # if kwargs has 'seed' key, use it, otherwise use default seed
    if kwargs.get('seed') is not None:
        pmpd.spec.randomSeed = kwargs.get('seed')
    # else:
    #     pmpd.spec.randomSeed = 3751 #to make the simulation fully reproducible
    if kwargs.get('MPI'): # if it is True, use MPI
        pmpd.mpiEnabled = True # This is essential as it enables the invocation of the MPI-parallelized ParaDRAM routines.

    pmpd.spec.overwriteRequested = False # overwrite old existing simulation files with the same name
    pmpd.spec.outputFileName = "./temp_sampler_out/" # the root-name of the output files
    pmpd.spec.chainSize = no_samples # the number of samples to be generated
    pmpd.spec.silentModeRequested = True
    if x_init is not None:
        pmpd.spec.startPointVec = x_init # the initial point of the Markov Chain
    pmpd.spec.targetAcceptanceRate = [0.1,0.3] # ensure the MCMC sampling efficiency does not become too large or too small.
    pmpd.spec.sampleSize = no_samples//5 # the final output sample size (optional)
    pmpd.runSampler ( ndim = n_dim, # number of dimensions``
                      getLogFunc = log_func   # the objective function
                    )
    pmpd.readSample() # generates pmpd.sampleList
    sample = pmpd.sampleList[0] # returns decorrelated samples. Its size is as in the *sample.txt file

    return sample.df
#%%
if __name__ == "__main__":
    # generate observed data
    # X = ss.norm(loc=3, scale=1).rvs(size=5000)

    # def guassian_posterior(theta):
    #     # returns the unnormalized log posterior
    #     loglik = np.sum(np.log(ss.norm(loc=theta, scale=1).pdf(X)))
    #     logprior = np.log(ss.norm(loc=0, scale=1).pdf(theta))
        
    #     return loglik + logprior


    # #%%
    # sample_df = MCMC_DRAM(guassian_posterior, n_dim=1)
    # sns.kdeplot((sample_df['SampleVariable1']))
    # plt.show()
    #%%
    XX = ss.multivariate_normal(mean=[3.0, 2.0], cov=[[0.5, 0.1], [0.1, 0.5]]).rvs(size=200)

    def MVN_posterior(theta, XX=XX, print_=False):
        cov_p = [[1.0, 0.5], [0.5, 1.0]]
       # cov_p = [[10.0, 8.5], [8.5, 10.0]]
        #cov_p = [[0.1, 0.01], [0.01, 0.1]]
        loglik = np.sum(np.log(ss.multivariate_normal(mean=theta, cov=[[0.5, 0.1], [0.1, 0.5]]).pdf(XX)))
        logprior = np.log(ss.multivariate_normal(mean=[1.0, 1.0], cov=cov_p).pdf(theta))

        return loglik + logprior
    def MVN_posterior_infer_log_noise(theta, XX=XX, print_=False):
        #cov_p = [[1.0, 0.5], [0.5, 1.0]]
        cov_p = [[10.0, 8.5], [8.5, 10.0]]
        #cov_p = [[0.1, 0.01], [0.01, 0.1]]
        loglik = np.sum(np.log(ss.multivariate_normal(mean=[theta[0],theta[1]], cov=[[np.exp(theta[2]), theta[3]], [theta[3], np.exp(theta[4])]]).pdf(XX)))
        logprior = np.log(ss.multivariate_normal(mean=[1.0, 1.0], cov=cov_p).pdf([theta[0],theta[1]]))

        return loglik + logprior
    
    # def tmp(theta, XX=XX):
    #     return MVN_posterior(XX, theta)
    #%%
    sample_df = MCMC_DRAM(MVN_posterior, n_dim=2, seed =666) #x_init=[3.0,2.0]
    #sample_df = MCMC_DRAM(MVN_posterior_infer_log_noise, n_dim=5, seed =666)
    #sample_df = MCMC_DRAM(tmp, n_dim=2)

    sns.kdeplot((sample_df['SampleVariable1']))
    plt.figure()
    sns.kdeplot((sample_df['SampleVariable2']))
    plt.show()

    # ------- different scale of the mean values
    #%%
    cov = [[1.0, 0.5], [0.5, 1.0]]
    # scale the above covariance matrix by 1E-02
    cov = np.array(cov) * 1E-02
    XX = ss.multivariate_normal(mean=[3.0, 2.0E-02], cov=cov).rvs(size=1)

    def MVN_posterior(theta, XX=XX, print_=False):
        loglik = np.sum(np.log(ss.multivariate_normal(mean=theta, cov=cov).pdf(XX)))
        logprior = np.log(ss.multivariate_normal(mean=[1.0, 1E-02], cov=[[1.0, 0.5], [0.5, 1.0]]).pdf(theta))

        return loglik + logprior
    # def tmp(theta, XX=XX):
    #     return MVN_posterior(XX, theta)
    #%%
    sample_df = MCMC_DRAM(MVN_posterior, n_dim=2)
    #sample_df = MCMC_DRAM(tmp, n_dim=2)

    sns.kdeplot((sample_df['SampleVariable1']))
    plt.figure()
    sns.kdeplot((sample_df['SampleVariable2']))
    plt.show()


# %%
