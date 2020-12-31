import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import yaml
from linear_model_error import LinearModelError
from bayes.multi_model_error import MultiModelError


class MultiLinearModelError(MultiModelError):
    """
    Define a class that combines multiple linear model errors
    """
    def __init__(self, yaml_file_list):
        """Create class combining different linear model errors

        Args:
            yaml_file_list(list of str): list of yaml files with the experimental data
        """
        super().__init__()
        for yaml_file in yaml_file_list:
            single_model_error = LinearModelError(str(yaml_file))
            # read in the parameters that are given in the experiment (in this case constant offset a)
            parameter = single_model_error.get_parameter_dict()
            # add parameters that are model parameters (not given in the experimental data)
            parameter.define("b")
            self.add(single_model_error, parameter)

        # define shared parameters (only one parameter b for all linear models)
        self.join(shared='b')
        # set this shared parameter to be latent (free parameters to be optimized)
        # no key given - shared variable
        self.set_latent('b')


class TestOptimize(unittest.TestCase):
    def test_linear_virtual_data(self):
        yaml_file_list_linear_experiment_data = \
            Path(Path(__file__).parents[0]).glob('virtual_linear_experiment_data_*.yaml')
        all_experiments_model_error = MultiLinearModelError(yaml_file_list_linear_experiment_data)

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
        all_experiments_model_error = MultiLinearModelError(yaml_file_list_quadratic_experiment_data)
        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        virtual_experiment_parameter_file = \
            Path(__file__).parent / 'virtual_quadratic_experiment_model.yaml'
        with open(virtual_experiment_parameter_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        # check to obtain linear coefficient using sympy
        # solve(diff(integrate(((c*x**2+b*x)-bbar*x)**2,(x,0,1)),bbar),bbar)
        # accuracy is limited due to integration error with a limited number of samples
        self.assertAlmostEqual(result.x[0], d['b']+3./4.*d['c'], 6)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
