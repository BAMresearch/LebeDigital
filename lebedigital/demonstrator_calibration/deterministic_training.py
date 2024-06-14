import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib as mpl
from matplotlib import rc

# set torch deafult data type to float32
torch.set_default_dtype(torch.float32)

#local imports
from lebedigital.demonstrator_calibration.parametric_model import NN_mean, train_NN

def indirect_train_NN(model:callable,x, z,latent_dim:int, forwardsolver:callable,
                       epochs=100, lr=1e-3, hidden_dim=20,forwardsolver_differentiable=True,**kwargs):
    """
    Parameters
    ----------
    model : callable
        model to be trained
    x : torch.tensor [N x D]
        input data (model input) N number of observed data and D dimensions
    z : torch.tensor [N x D]
        output data (physcis based model output) N number of observed data and D dimensions
    latent_dim : int
        latent dimension of the model (dim(b))
    forwardsolver : callable
        callable that outputs of the model and returns y. Needs kwargs to be passed.
    epochs : int, optional
        number of epochs, by default 1000
    lr : float, optional
        learning rate, by default 1e-3
    hidden_dim : int, optional
        hidden dimension of the NN, by default 512

    Returns
    -------
    nn_mean : NN_mean
        trained NN_mean model
    """
    input_dim = x.shape[1]
    #output_dim = z.shape[1]
    output_dim = latent_dim
    # define the model
    nn_mean = model(input_dim, hidden_dim, output_dim)
    # define the loss function
    criterion = nn.MSELoss()
    # define the optimizer
    optimizer = optim.Adam(nn_mean.parameters(), lr=lr)
    
    # train the model
    for epoch in range(epochs):
        # forward pass of the NN
        b_pred = nn_mean(x)

        if forwardsolver_differentiable:
            # foward solver as an observation operator
            z_pred = forwardsolver(b_pred,**kwargs)

            # compute the loss
            loss = criterion(z_pred, z)
        else:
            # take b_pred as mean of a MV Gaussian distribution in torch with a diagonal covariance matrix
            b_pred_dist = torch.distributions.MultivariateNormal(b_pred, torch.eye(latent_dim))

            # define a log likelihood function using MV gaussian distribution in torch
            def log_likelihood(b):
                
            




        # backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # save the NN gradient norm for all the wieghts and biases
        nn_grad_norm = torch.norm(torch.cat([p.grad.flatten() for p in nn_mean.parameters()]))

        if epoch % 100 == 0:
            print(f"Epoch {epoch}, Loss {loss.item():.4f}, NN grad norm {nn_grad_norm.item():.4f}")
    return nn_mean

# function to do prediction for the validation data with the trained NN
def predict_NN(nn_mean:callable, x, z, forwardsolver:callable, **kwargs):
    """
    Parameters
    ----------
    nn_mean : callable
        trained NN_mean model
    x : torch.tensor [N x D]
        input data (model input) N number of observed data and D dimensions
    z : torch.tensor [N x D]
        output data (physcis based model output) N number of observed data and D dimensions
    forwardsolver : callable
        callable that outputs of the model and returns y. Needs kwargs to be passed.

    Returns
    -------
    y_pred : torch.tensor [N x D]
        predicted output data (model output) N number of observed data and D dimensions
    """
    # forward pass of the NN
    b_pred = nn_mean(x)

    # foward solver as an observation operator
    z_pred = forwardsolver(b_pred,**kwargs)

    # prediction accuracy
    loss = torch.norm(z_pred - z)

    print(f"Prediction accuracy: {loss.item():.4f}")
    return z_pred, loss