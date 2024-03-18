import pytest
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
import matplotlib as mpl
from matplotlib import rc

from lebedigital.demonstrator_calibration.parametric_model import NN_mean, train_NN
# set torch deafult data type to float32
torch.set_default_dtype(torch.float32)
# set seed for reproducibility
torch.manual_seed(0)

datetime = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")

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

def test_train_NN():
    x = torch.tensor([[0.3],[0.6]])
    #y = torch.tensor([[2.916E-4, 0.0024229, 5.554, 500e3]])
    y = torch.tensor([[2.916, 2.4229, 5.554, 5.0],[2.7, 2.43, 5.56, 4.8]])
    nn_mean = train_NN(NN_mean,x, y, epochs=400, lr=1e-2, hidden_dim=10)
    assert nn_mean(x).shape == (2, 4)

    x_test = torch.tensor([[0.1], [0.2]])
    y_pred = nn_mean(x_test)
    assert y_pred.shape == (2, 4)
    y_true = torch.tensor([[2.7957, 2.4180, 5.5067, 4.8626],
        [2.8008, 2.4217, 5.5266, 4.8775]])
    # assert the y_pred and y_true are close
    print(f'y_pred = {y_pred}')
    assert torch.allclose(y_pred, y_true, rtol=1e-3, atol=1e-3)

#test_NN_mean()
#test_train_NN()