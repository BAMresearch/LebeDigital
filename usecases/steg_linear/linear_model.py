import numpy as np


class LinearModel:
    """
    Define a standard forward model for a specific type of test (in this case a linear model)

    Note:
        f(x) = a + b*x , where a is a measured deterministic parameter, b is a model parameter
        to be calibrated and x are the sensor positions in the interval [0,1].
        Different sensor positions for the function f(x) as well as for the derivative
        f'(x) = b can be prescribed.

    Attributes:
        x_function (np.array): Array of positions x where function values are measured
        x_derivative (np.array): Array of positions x where function derivatives are measured
    """

    def __init__(self, x_function, x_derivative):
        """Construct a linear forward model with sensor positions for function values and derivatives"""
        self.x_function = x_function
        self.x_derivative = x_derivative

    def __call__(self, parameters):
        """Evaluate the linear model at predefined sensor positions for function values and derivatives
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (a and b)

        Returns:
            [f_x, df_x]: function values and derivatives at the locations defined in the constructor
        """
        self.check_parameters(parameters)
        f_x = parameters['a'] + parameters['b'] * self.x_function
        df_x = parameters['b'] * np.ones(len(self.x_derivative))
        return [f_x, df_x]

    @staticmethod
    def are_parameters_valid(parameters):
        """
        Check the parameters (all required are given and are also within the bounds)
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes
        Returns:
            true: if check is successful, other wise return false
        """
        try:
            LinearModel.check_parameters()
        except Exception as e:
            print(e)
            return False
        return True

    @staticmethod
    def check_parameters(parameters):
        """
        Check the parameters (all required are given and are also within the bounds)
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes
        Raises:
            throws an exception when parameters are either not given or
            out of bounds
        """
        assert 'a' in parameters.names
        assert 'b' in parameters.names

        assert parameters['a'] > 0
        assert parameters['b'] > 0
