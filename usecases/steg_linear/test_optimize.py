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
    def __init__(self, yaml_file_metadata, yaml_file_data):
        """Create class combining different linear model errors

        Note:
            Currently implemented as a class to later add individual routines of adding additional errors
            e.g. using database queries, etc. currently, everything is only in yaml files and thus in the
            constructor.

        Args:
            yaml_file_list(list of str): list of yaml files with the experimental data
        """
        super().__init__()

        with open(yaml_file_metadata, "r") as f:
            meta_data = yaml.load(f, Loader=yaml.FullLoader)

        x_function = np.asarray(meta_data['x_function'])
        x_derivative = np.asarray(meta_data['x_derivative'])
        all_a = np.asarray(meta_data['all_a'])

        with open(yaml_file_data, "r") as f:
            data = yaml.load(f, Loader=yaml.FullLoader)

        for idx, a in enumerate(all_a):
            data_f = np.asarray(data[idx]['f'])
            data_df = np.asarray(data[idx]['df'])
            assert(a == data[idx]['a'])
            single_model_error = LinearModelError(x_function, x_derivative, data_f, data_df, a)
            # read in the parameters that are given in the experiment (in this case constant offset a)
            parameter = single_model_error.get_parameter_dict()
            # add parameters that are model parameters (not given in the experimental data)
            parameter.define("b")
            self.add(single_model_error, parameter)

        # define one shared latent parameter (only one parameter b for all linear models)
        self.latent.add_by_name('b')


class TestOptimize(unittest.TestCase):
    def test_linear_virtual_data(self):   # data without noise
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_model_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_linear_model_meta.yaml',
                                  'virtual_experiment_linear_model_data.yaml')

        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        #print('steg, results, TestOptimize, least_squres:',result) # steg

        # assert quadratic component of virtual experiment to be zero
        self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        self.assertAlmostEqual(result.x[0], d['b'], 7) # 7? steg
        # check cost function to be zero
        self.assertAlmostEqual(result.cost, 0)

    def test_linear_virtual_data_split(self):   # data with noise
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_model_with_noise_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_linear_model_with_noise_meta.yaml',
                                  'virtual_experiment_linear_model_with_noise_data.yaml')

        def weighted_all_experiments_model_error(prm_vector):
            all_model_errors = all_experiments_model_error.evaluate(prm_vector)
            obj = 0.
            for key, me in all_model_errors.items():
                obj = obj + np.dot(me['f'], me['f']) / d['sigma_noise_function'] ** 2 + \
                      np.dot(me['df'], me['df']) / d['sigma_noise_derivative'] ** 2
            return np.sqrt(obj)

        weighted_start_vector = np.array([0.7])
        weighted_result = least_squares(weighted_all_experiments_model_error, weighted_start_vector)

        # assert quadratic component of virtual experiment to be zero
        self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        self.assertAlmostEqual(weighted_result.x[0], d['b'], 4) # 4? steg 
'''
    def test_quadratic_virtual_data(self):
        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_quadratic_model_meta.yaml',
                                  'virtual_experiment_quadratic_model_data.yaml')
        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_quadratic_model_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        # check to obtain linear coefficient using sympy
        # solve(diff(integrate(((c*x**2+b*x)-bbar*x)**2,(x,0,1)),bbar),bbar)
        # accuracy is limited due to integration error with a limited number of samples
        self.assertAlmostEqual(result.x[0], d['b']+3./4.*d['c'], 6)
'''

if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
