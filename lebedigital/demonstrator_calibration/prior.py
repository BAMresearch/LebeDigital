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

# setting seed for reproducibility in torch
th.manual_seed(0)

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
        In the subsequent, x is the input parameters and b are the latent parameters with (m x d), 
        m being the number of samples and d being the dimension of the latent space.

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
        # check if self.cov_params is a tensor and requires grad
        if not isinstance(self.cov_params, th.Tensor):
            self.para_cov_torch = th.tensor(self.cov_params,requires_grad=True)
        else:
            self.para_cov_torch = self.cov_params
        
    def cov(self):
        # compute the covariance matrix
        # the parameter are positioned as =(0,0), (1,0), (1,1), (2,0), (2,1), (2,3) ...
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
            #L[np.diag_indices(self.latent_dim)] = th.exp(L[np.diag_indices(self.latent_dim)])
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
        """_summary_

        Parameters
        ----------
        x : _type_
            _description_
        b : list with size (1, latent_dim)
            _description_

        Returns
        -------
        _type_
            _description_
        """
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
        if x.requires_grad:
            grad_mean , grad_cov = x.grad, self.para_cov_torch.grad
        else:
            # get grad of the nn which is the self.mean parameters
            #grad_mean = th.cat([p.grad.flatten() for p in self.mean.parameters()])
            grad_mean = [p.grad.detach().clone() for p in self.mean.parameters()]
            grad_cov = self.para_cov_torch.grad.detach().clone()
        # return the detched and cloned gradeints
        #return grad_mean.detach().clone(), grad_cov.detach().clone()
        return grad_mean, grad_cov
    
    def grad_estimate_score_function(self, x:np.array,b:list, return_grad:bool=False): 
        """Grad estimate using the score function estimator for the M-step of the EM algorithm.

        Parameters
        ----------
        x : array
            array with N x input_dim, with N being number of observed data points
        b : list
            list of len N of arrays of size (n_samples, latent_dim)
        return_grad : bool, optional
            flag to say return the grad itself or the expected log_prob to be later used for computing the grads, by default False
        """
        # compute the gradient of the log pdf
        if not isinstance(x, th.Tensor):
            x = th.tensor(x)
        if not isinstance(b, th.Tensor):
            #b = th.tensor(b)
            b = [th.tensor(b[i]) for i in range(len(b))]
        assert len(x.shape) == 2, "x should be a 2d tensor"
        assert x.shape[0] == len(b), "The number of x and b should be the same"
        grad_mean, grad_cov, log_prob = [], [], []
        for i in range(len(b)): # iterating over the data pairs
            for m in range(b[i].shape[0]): # iterating over the numbvber of samples
                if return_grad:
                    grad_mean_, grad_cov_ = self.grad_log_pdf(x[i,:],b[i][m,:])
                    grad_mean.append(grad_mean_)
                    grad_cov.append(grad_cov_)
                else:
                    # compute mean of the log_prob
                    log_prob.append(self.log_prob(x[i,:],b[i][m,:]))

        if return_grad:
            #expected_grad_mean = th.mean(th.stack(grad_mean),dim=0)
            expected_grad_mean = self._mean_nn_grad(grad_mean)
            expected_grad_cov = th.mean(th.stack(grad_cov),dim=0)
        else:
            expected_log_prob = th.mean(th.stack(log_prob),dim=0)
    
        # return the grad itself or the expected log_prob
        if return_grad:
            return expected_grad_mean, expected_grad_cov
        else:
            return expected_log_prob

    def plot(self,x, n_samples):
        samples = self.sample(x,n_samples)
        #sb.kdeplot(samples[:,0], samples[:,1], shade=True, cmap='Blues')
        plt.plot(samples[:,0], samples[:,1], 'o')
        plt.show()
    
    def _mean_nn_grad(self,list_:list)->list:
        """_summary_

        Parameters
        ----------
        list_ : list of lists with each list i being the NN parameter grad for the ith data point, j being the weights for a biases for a respective layer
            _description_

        Returns
        -------
        list
            mean of the grad of the NN parameters
        """
        grad_mean = [th.mean(th.stack([list_[j][i] for j in range(len(list_))]),dim=0) for i in range(len(list_[0]))]

        return grad_mean



if __name__ == "__main__":

# ----------------------------------------------------------------------------------------
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
        nn_mean = train_NN(NN_mean,x, y, epochs=2000, lr=1e-2, hidden_dim=10)

        cov_params = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

        prior_ = prior(nn_mean, cov_params=cov_params, cov_type='full',latent_dim=4)
        cov = prior_.cov()
        print(f'the covariance matrix is {cov}')
        log_prob = prior_.log_prob([0.4],[2.7, 2.43, 5.56, 4.8])
        print(f'the log prob is {log_prob}')
        sample = prior_.sample([0.4],1000)
        print(f'the sample mean is {np.mean(sample,axis=0)}')
        g_mean, g_cov = prior_.grad_log_pdf([0.6],[2.7, 2.43, 5.56, 4.8])
        print(f'the gradient of mean is {g_mean} and the gradient of cov is {g_cov}')

    test_prior_with_nn()

    #%%
    # run the pretraining for 1 dim input and 4 dim output synthetic data 
    x = th.tensor([[0.3],[0.6]])
    #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
    y = th.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
    nn_mean = train_NN(NN_mean,x, y, epochs=2000, lr=1e-2, hidden_dim=10)

    cov_params = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

    prior_ = prior(nn_mean, cov_params=cov_params, cov_type='full',latent_dim=4)
    cov = prior_.cov()
    print(f'the covariance matrix is {cov}')
    #log_prob = prior_.log_prob([0.4],[np.array([[2.7, 2.43, 5.56, 4.8],[2.72, 2.41, 5.50, 4.84]])])
    #g_mean, g_cov = prior_.grad_log_pdf([0.4],[np.array([[2.7, 2.43, 5.56, 4.8],[2.72, 2.41, 5.50, 4.84]])])

    # create a list of len 4 with 2d array of size 2x4
    b_list = [np.array([[2.71, 2.45, 5.5, 4.84],[2.75, 2.2, 5.1, 4.80]]),
            np.array([[2.7, 2.43, 5.56, 4.8],[2.72, 2.41, 5.50, 4.84]]),
            np.array([[2.6, 2.33, 5.26, 4.5],[2.22, 2.41, 5.10, 4.24]]),
                np.array([[2.1, 2.13, 5.26, 4.4],[2.12, 2.21, 5.40, 4.54]])]

    # expected_log_prob = prior_.grad_estimate_score_function([[0.2],[0.35],[0.4],[0.5]],b_list)
    # expected_log_prob.backward()
    # grad_mean = th.cat([p.grad.flatten() for p in prior_.mean.parameters()])
    # grad_cov = prior_.para_cov_torch.grad

    # grad_direct_mean, grad_direct_cov = prior_.grad_estimate_score_function([[0.2],[0.35],[0.4],[0.5]],b_list,return_grad=True)

    # values at which the gradient is expected to be almost close to zero
    b_list_new = [np.array([[2.7, 2.43, 5.56, 4.8],[2.7, 2.43, 5.56, 4.8],[2.7, 2.43, 5.56, 4.8]])]
    grad_direct_mean, grad_direct_cov = prior_.grad_estimate_score_function([[0.6]],b_list_new,return_grad=True)


    #%%
    temp_cov_param = th.tensor([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],requires_grad=True)
    # pass nn_mean.parameters() and temp_cov_param to the optimizer
    optimizer = th.optim.Adam([{'params': nn_mean.parameters()},{'params': temp_cov_param}], lr=1e-2)
    for i,p in enumerate(nn_mean.parameters()):
        p.grad = g_mean[i]
    optimizer.step()
    for i,p in enumerate(nn_mean.parameters()):
        print(p)




