import numpy as np
from linear_model import LinearModel
from bayes.parameters import ModelParameters
from collections import OrderedDict


class LinearModelError:
    """
    Define a standard model error for a specific type of test (in this case a linear model)

    Note:
        f(x) = a + b*x , where a is a measured deterministic parameter, b is a model parameter
        to be calibrated and x are the sensor positions in the interval [0,1].
        Different sensor positions for the function f(x) as well as for the derivative
        f'(x) = b can be prescribed.

    Attributes:
        x_function (np.array): Array of positions x where function values are measured
        x_derivative (np.array): Array of positions x where function derivatives are measured
        data_f (np.array): Array of sensor readings at x_function
        data_df (np.array): Array of sensor reading at x_derivative
        a (float): offset of the linear function (scalar deterministic value)
        linear_model (LinearModel): linear forward model (a+b*x) , where 'b' is a model parameter
    """

    def __init__(self, x_function, x_derivative, data_f, data_df, a):
        """Create a linear model error

        Args:
            x_function (np.array): Array of positions x where function values are measured
            x_derivative (np.array): Array of positions x where function derivatives are measured
            data_f (np.array): Array of sensor readings at x_function
            data_df (np.array): Array of sensor reading at x_derivative
            a (float): offset of the linear function (scalar deterministic value)
        """
        self.x_function = x_function
        self.x_derivative = x_derivative
        self.data_f = data_f
        self.data_df = data_df
        self.a = a
        self.linear_model = LinearModel(self.x_function, self.x_derivative)

    def __call__(self, parameters):
        """Evaluate the model error for a specific experiment

        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (a and b)

        Returns:
            concatenated model error (np.array) of all individual model errors
        """
        [f, df] = self.linear_model(parameters)
        return np.concatenate((f-self.data_f, df-self.data_df))

    def evaluate(self, parameters):
        """Evaluate the model error for a specific experiment

        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (a and b)

        Returns:
            return an OrderDict of the individual model errors (e.g. per sensor)
        """
        [f, df] = self.linear_model(parameters)
        return OrderedDict([('f', f-self.data_f), ('df', df-self.data_df)])

    def get_parameter_dict(self):
        """Create a parameter list initialized with parameters given in the experimental data file

        Returns: parameter list
        """
        prm = ModelParameters()
        prm.define("a", self.a)
        return prm
