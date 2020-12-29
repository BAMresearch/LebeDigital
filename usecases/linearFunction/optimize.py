import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import yaml
from linear_model_error import *
from bayes.parameters import *
from bayes.multi_model_error import MultiModelError


class AllExperimentsModelError(MultiModelError):
    "concatenate multiple 'standardized' model errors into a joint class and take care of the shared parameters"

    def __init__(self, yaml_file_list):
        # For the inference, we combine them and use a 'key' to distinguish
        # e.g. "A" from the one model to "A" from the other one.
        super().__init__()
        for yaml_file in yaml_file_list:
            single_model_error = LinearModelError(str(yaml_file))
            parameter = single_model_error.get_parameter_dict()
            parameter.define("b")
            key = self.add(single_model_error, parameter)

        self.join(shared='b')
        self.set_latent('b')


class TestOptimize(unittest.TestCase):
    def test_linear_virtual_data(self):
        yaml_file_list_linear_experiment_data = \
            Path(Path(__file__).parents[0]).glob('virtual_linear_experiment_data_*.yaml')
        all_experiments_model_error = AllExperimentsModelError(yaml_file_list_linear_experiment_data)

        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        virtual_experiment_parameter_file = \
            Path(__file__).parent / 'virtual_linear_experiment_model.yaml'
        with open(virtual_experiment_parameter_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        # assert quadratic component of virtual experiment to be zero
        self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        self.assertAlmostEqual(result.x[0], d['b'])
        # check cost function to be zero
        self.assertAlmostEqual(result.cost, 0)

    def test_quadratic_virtual_data(self):
        yaml_file_list_quadratic_experiment_data = \
            Path(Path(__file__).parents[0]).glob('virtual_quadratic_experiment_data_*.yaml')
        all_experiments_model_error = AllExperimentsModelError(yaml_file_list_quadratic_experiment_data)
        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        #print(f"optimal parameter {result.x} with cost function {result.cost}.")

        virtual_experiment_parameter_file = \
            Path(__file__).parent / 'virtual_quadratic_experiment_model.yaml'
        with open(virtual_experiment_parameter_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        # check to obtain linear coefficient using sympy
        # solve(diff(integrate(((c*x**2+b*x)-bbar*x)**2,(x,0,1)),bbar),bbar)
        # accuracy is limited due to integration error with limited number of samples
        self.assertAlmostEqual(result.x[0], d['b']+3./4.*d['c'], 6)

if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
