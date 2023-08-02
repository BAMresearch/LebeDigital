#%%
import torch as th
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

from datetime import datetime
import matplotlib as mpl
from matplotlib import rc

from lebedigital.demonstrator_calibration.parametric_model import NN_mean, train_NN

# set torch deafult data type to float32
th.set_default_dtype(th.float64)

class prior:
    def __init__(self, mean:callable, cov_params:list, latent_dim:int, cov_type:str='diag'):
        """        This class defines the prior distribution of the parameters
        of the forward model. The prior is a multivariate normal distribution
        with mean and covariance matrix. The mean is a function of the input
        parameters (NN for eg) and the covariance matrix is a function of the covariance
        parameters. The covariance matrix can be either diagonal or full.
        The covariance parameters are the parameters of the covariance matrix.
        The covariance matrix is computed as follows:
        1. If the covariance matrix is diagonal, then the covariance parameters
        are the diagonal elements of the covariance matrix. The diagonal elements
        are exponentiated to ensure positivity.
        2. If the covariance matrix is full, then the covariance parameters are
        the lower triangular elements of the covariance matrix. The lower triangular
        elements are exponentiated to ensure positivity. The diagonal elements are
        exponentiated to ensure positivity and then the lower triangular matrix is
        computed using the Cholesky decomposition.

        Parameters
        ----------
        mean : instance of nn.Module
            Mostly its the NN, not implemebted for other cases yet.
        cov_params : list
            parameters of the covariance matrix. If the covariance matrix is diagonal,
            then the covariance parameters are the diagonal elements of the covariance
            matrix. The diagonal elements are exponentiated to ensure positivity.
            If the covariance matrix is full, then the covariance parameters are
            the lower triangular elements of the covariance matrix. 
        latent_dim : int
            The latents which are the parameters of the forward model.
        cov_type : str, optional
            accepts 'full' or 'diag' by default 'diag'
        """



        self.mean = mean # here goes the NN
        self.cov_params = cov_params
        self.cov_type = cov_type
        self.latent_dim = latent_dim
        self.para_cov_torch = th.tensor(self.cov_params,requires_grad=True)
        
    def cov(self):
        # compute the covariance matrix
        #para = th.tensor(self.cov_params,requires_grad=True)
        if self.cov_type == 'diag':
            return th.diag(th.exp(self.para_cov_torch))
        elif self.cov_type == 'full':
            # for N dim matrix, check the number of elements in the lower triangular matrix
            # if it is N*(N+1)/2 then it is a lower triangular matrix
            assert len(self.cov_params) == self.latent_dim*(self.latent_dim+1)/2,\
            "The number of parameters for the covariance matrix is not correct"
            # compute the lower triangular matrix
            L = th.zeros(self.latent_dim,self.latent_dim)
            L[np.tril_indices(self.latent_dim)] = self.para_cov_torch
            # diagonal elements are exponentiated
            L[np.diag_indices(self.latent_dim)] = th.exp(L[np.diag_indices(self.latent_dim)])
            # return the covariance matrix
            return th.mm(L,L.t())
        else:
            raise ValueError("The covariance type is not correct. It should be either diag or full")

    def sample(self,x, n_samples):
        # convert x to tensor if it is not and it should be atleast 1d
        #x = np.atleast_1d(x)
        if not isinstance(x, th.Tensor):
            x = th.tensor(x)
        return th.distributions.MultivariateNormal(self.mean(x), self.cov()).sample((n_samples,)).detach().numpy()

    def log_prob(self,x, b):
        # convert x and b to tensor if it is not
        #x = np.atleast_1d(x)
        #b = np.atleast_1d(b)
        if not isinstance(x, th.Tensor):
            x = th.tensor(x)
        if not isinstance(b, th.Tensor):
            b = th.tensor(b)
        return th.distributions.MultivariateNormal(self.mean(x), self.cov()).log_prob(b)
    
    def grad_log_pdf(self,x, b):
        # convert x and b to tensor if it is not
        #x = np.atleast_1d(x)
        #b = np.atleast_1d(b)
        if not isinstance(x, th.Tensor):
            x = th.tensor(x)
        if not isinstance(b, th.Tensor):
            b = th.tensor(b)
        # compute the gradient of the log pdf
        log_pdf = th.distributions.MultivariateNormal(self.mean(x), self.cov()).log_prob(b)
        log_pdf.backward()
        # return the gradient of mean and covariance w.r.t the parameters
        # TODO: should be grad wrt NN parameters
        if x.requires_grad:
            grad_mean , grad_cov = x.grad, self.para_cov_torch.grad
        else:
            # get grad of the nn which is the self.mean parameters
            grad_mean = th.cat([p.grad.flatten() for p in self.mean.parameters()])
            grad_cov = self.para_cov_torch.grad

        return np.array(grad_mean), np.array(grad_cov)


    def plot(self,x, n_samples):
        samples = self.sample(x,n_samples)
        #sb.kdeplot(samples[:,0], samples[:,1], shade=True, cmap='Blues')
        plt.plot(samples[:,0], samples[:,1], 'o')
        plt.show()

#%%
# writre a test for all the class methods above with mean being a 2*x function
def test_prior():
    def mean(x):
        return 2*x
    # define x and cov_params
    x = th.tensor([1.0,1.0],requires_grad=True)
    cov_params = [0.01,0.01,0.01]

    prior_ = prior(mean, cov_params=cov_params, cov_type='full',latent_dim=2)
    cov = prior_.cov()
    print(f'the covariance matrix is {cov}')
    sample = prior_.sample(x,1000)
    g_mean, g_cov = prior_.grad_log_pdf(x, [2.0,2.0])
    log_prob = prior_.log_prob([1.0,1.0],[2.0,2.0])
    print(f'the gradient of mean is {g_mean} and the gradient of cov is {g_cov}')
    print(f'the sample mean is {np.mean(sample,axis=0)}and the log prob is {log_prob}')
    #prior_.plot([1.0,1.0],1000)

#%%
#test_prior()

#%%
def test_prior_with_nn():
    # run the pretraining for 1 dim input and 4 dim output synthetic data 
    x = th.tensor([[0.3],[0.6]])
    #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
    y = th.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
    nn_mean = train_NN(NN_mean,x, y, epochs=800, lr=1e-2, hidden_dim=10)

    cov_params = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

    prior_ = prior(nn_mean, cov_params=cov_params, cov_type='full',latent_dim=4)
    cov = prior_.cov()
    print(f'the covariance matrix is {cov}')
    log_prob = prior_.log_prob([0.4],[2.7, 2.43, 5.56, 4.8])
    print(f'the log prob is {log_prob}')
    sample = prior_.sample([0.4],1000)
    print(f'the sample mean is {np.mean(sample,axis=0)}')
    g_mean, g_cov = prior_.grad_log_pdf([0.4],[2.7, 2.43, 5.56, 4.8])
    print(f'the gradient of mean is {g_mean} and the gradient of cov is {g_cov}')

test_prior_with_nn()
