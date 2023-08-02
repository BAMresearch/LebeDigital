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

datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

# class for NN based mean 
# TODO : inherit it from base class of paramteirc models. Can be as simple as 
# linear regression also
class NN_mean(nn.Module):
    # TODO: add option to specity depth of the NN
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(NN_mean, self).__init__()
        self.layer1 = nn.Linear(input_dim, hidden_dim)
        self.layer2 = nn.Linear(hidden_dim, hidden_dim)
        #self.layer3 = nn.Linear(512, 512)
        self.layer3 = nn.Linear(hidden_dim, output_dim)
        #self.dropout = nn.Dropout(p=0.2)

    def forward(self, x):
        """_summary_

        Parameters
        ----------
        x : tensor
            Input feature vector with N data points and D dimensions

        Returns
        -------
        b : tensor [NxD]
            Output feature vector with N data points and D dimensions
            _description_
        """
        x = torch.relu(self.layer1(x))
        #x = self.dropout(x)
        x = torch.relu(self.layer2(x))
        #x = self.dropout(x)
        #x = torch.relu(self.layer3(x))
        x= self.layer3(x)
        #x = self.dropout(x)
        #x = self.layer4(x)
        return x

# write a pytest for the above
def test_NN_mean():
    """_summary_
    """
    # create a dummy input
    input_dim = 10
    hidden_dim = 20
    output_dim = 5
    x = torch.rand(100, input_dim)
    # create a dummy model
    model = NN_mean(input_dim, hidden_dim, output_dim)
    # check the output size
    assert model(x).shape == (100, output_dim)

#function to overload the parameters of the NN_mean by a prescribed value
def overload_params(model, params):
    """_summary_

    Parameters
    ----------
    model : nn.Module
        pytorch model
    params : list
        list of parameters to overload

    Returns
    -------
    model : nn.Module
        pytorch model with overloaded parameters
    """
    # get the state dictionary of the model
    state_dict = model.state_dict()
    # loop over the parameters to overload
    for key, value in params.items():
        # overload the parameter
        state_dict[key] = value
    # load the state dictionary back to the model
    model.load_state_dict(state_dict)
    return model


# pretrain the above model to get a good initialization
def train_NN(model:callable,x, y, epochs=100, lr=1e-3, hidden_dim=20):
    """
    Parameters
    ----------
    model : callable
        model to be trained
    x : torch.tensor [N x D]
        input data N number of observed data and D dimensions
    y : torch.tensor [N x D]
        output data N number of observed data and D dimensions
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
    output_dim = y.shape[1]

    # define the model
    # check if the model neneds to be initialized
    if isinstance(model, nn.Module):
        nn_mean = model   # if pre trained model is passed
    else:
        nn_mean = model(input_dim, hidden_dim, output_dim)
    # define the loss function
    criterion = nn.MSELoss()

    # define the optimizer
    optimizer = optim.Adam(nn_mean.parameters(), lr=lr)

    # train the model
    for epoch in range(epochs):
        # forward pass
        y_pred = nn_mean(x)
        loss = criterion(y_pred, y)

        # backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # save the NN gradient norm for all the wieghts and biases
        nn_grad_norm = torch.norm(torch.cat([p.grad.flatten() 
                                             for p in nn_mean.parameters()]))

        if epoch % 100 == 0:
            print(f"Epoch {epoch}, Loss {loss.item():.4f}, NN grad norm {nn_grad_norm.item():.4f}")
    
    # print the forward pass
    #y_pred = nn_mean(x)
    #print(f'predicted output: {y_pred}')
    return nn_mean
# write aa pytest for the above
def test_train_NN():
    """_summary_
    """
    # create a dummy input
    input_dim = 10
    hidden_dim = 20
    output_dim = 5
    x = torch.rand(100, input_dim)
    y = torch.rand(100, output_dim)
    # create a dummy model
    model = NN_mean(input_dim, hidden_dim, output_dim)
    # train the model
    model = train_NN(model, x, y, epochs=100, lr=1e-3, hidden_dim=20)
    # check the output size
    assert model(x).shape == (100, output_dim)
if __name__=='__main__':
# ------------ pre training ---------------------

    # run the pretraining for 1 dim input and 4 dim output synthetic data 
    x = torch.tensor([[0.3],[0.6]])
    #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
    y = torch.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
    nn_mean = train_NN(NN_mean,x, y, epochs=400, lr=1e-2, hidden_dim=10)

    #nn_mean = train_NN(nn_mean_1,x, y, epochs=400, lr=1e-2, hidden_dim=10)
    # predict for 4 different values of x irnage 0.1 to 0.8 with the trained model
    x_test = torch.tensor([[0.1], [0.2], [0.3], [0.4], [0.5], [0.6]])
    y_pred = nn_mean(x_test)
    print(f'predicted output: {y_pred}')



# test overload parameters function with the paramters of the pre-trained model




