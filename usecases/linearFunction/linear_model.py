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
        f_x = parameters['a'] * np.ones(len(self.x_function)) + parameters['b'] * self.x_function
        df_x = parameters['b'] * np.ones(len(self.x_derivative))
        return [f_x, df_x]

    @staticmethod
    def check_parameters(parameters, runtime_error=True):
        """
        Check the parameters (all required are given and are also within the bounds)
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes
            runtime_error(bool): specify if a RuntimeError is thrown (raise=True) or only
                         false returned on failed parameter checks
        Raises:
            throws a RuntimeError if raise=True when parameters are either not given or
            out of bounds
        Returns:
            true: if check is successful, other wise return false (or throw)
        """
        if 'a' not in parameters.names:
            if runtime_error is True:
                raise RuntimeError(
                    "Parameter 'a' not given as input to LinearModel."
                )
            else:
                return False
        if (parameters['a'] <= 0):
            if runtime_error is True:
                raise RuntimeError(
                    f"Parameter 'a' in LinearModel must be positive (={parameters['a']}"
                )
            else:
                return false

        if 'b' not in parameters.names:
            if runtime_error is True:
                raise RuntimeError(
                    "Parameter 'b' not given as input to LinearModel."
                )
            else:
                return False
        if (parameters['b'] <= 0):
            if runtime_error is True:
                raise RuntimeError(
                    f"Parameter 'b' in LinearModel must be positive (={parameters['b']}"
                )
            else:
                return False
        return True
