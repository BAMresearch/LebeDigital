import numpy as np
from linear_model import LinearModel
from bayes.parameters import *


class LinearModelError:
    def __init__(self, experimental_data_file):
        import yaml

        with open(experimental_data_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        self.x_function = np.asarray(d['x_function'])
        self.x_derivative = np.asarray(d['x_derivative'])
        self.data_f = np.asarray(d['f'])
        self.data_df = np.asarray(d['df'])
        self.a = d['a']
        self.linear_model = LinearModel(self.x_function, self.x_derivative)

    def __call__(self, parameters):
        [f, df] = self.linear_model(parameters)
        return np.concatenate((f-self.data_f, df-self.data_df))

    def get_parameter_dict(self):
        prm = ModelParameters()
        prm.define("a", self.a)
        return prm





