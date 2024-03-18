import pytest
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt

from datetime import datetime
import matplotlib as mpl
from matplotlib import rc

from lebedigital.demonstrator_calibration.sampler import MCMC_DRAM
# set torch deafult data type to float32
torch.set_default_dtype(torch.float32)
# set seed for reproducibility
torch.manual_seed(0)

datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

# designed to work locally, dont know how to install paramonte remotely. 
# install notes : https://www.cdslab.org/paramonte/generic/latest/installation/QUICKSTART.md
def test_MCMC_DRAM():
    # generate observed data
    XX = ss.multivariate_normal(mean=[3.0, 2.0], cov=[[0.5, 0.1], [0.1, 0.5]]).rvs(size=200)

    # define lkl and prior
    def MVN_posterior(theta, XX=XX, print_=False):
        cov_p = [[1.0, 0.5], [0.5, 1.0]]
       # cov_p = [[10.0, 8.5], [8.5, 10.0]]
        #cov_p = [[0.1, 0.01], [0.01, 0.1]]
        loglik = np.sum(np.log(ss.multivariate_normal(mean=theta, cov=[[0.5, 0.1], [0.1, 0.5]]).pdf(XX)))
        logprior = np.log(ss.multivariate_normal(mean=[1.0, 1.0], cov=cov_p).pdf(theta))
        return loglik + logprior

    sample_df = MCMC_DRAM(MVN_posterior, n_dim=2, seed =666, x_init=[3.0,2.0])

    # assert the mean
    assert np.mean(sample_df['SampleVariable1']) == pytest.approx(3.0, abs=0.1)
    assert np.mean(sample_df['SampleVariable2']) == pytest.approx(2.0, abs=0.1)

