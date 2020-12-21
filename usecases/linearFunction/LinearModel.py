import numpy as np

class LinearModel:
    def __init__(self, x_function, x_derivative):
        self.x_function = x_function
        self.x_derivative = x_derivative

    def  __call__(self, parameters):
        self.check_parameters(parameters)
        f_x =  parameters['a'] * np.ones(len(self.x_function)) + \
               parameters['b'] * self.x_function
        df_x = parameters['b'] * np.ones(len(self.x_derivative))
        return [f_x, df_x]

    def check_parameters(self, parameters):
        if 'a' not in parameters.names:
            raise RuntimeError(
                "Parameter 'a' not given as input to LinearModel."
            )
            if (parameters['a'] <= 0):
                raise RuntimeError(
                    f"Parameter 'a' in LinearModel must be positive (={parameters['a']}"
                )

        if 'b' not in parameters.names:
            raise RuntimeError(
                "Parameter 'b' not given as input to LinearModel."
            )
            if (parameters['b'] <= 0):
                raise RuntimeError(
                    f"Parameter 'b' in LinearModel must be positive (={parameters['b']}"
                )





