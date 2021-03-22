import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import yaml
from asymptotic_model_error import AsymptoticModelError
from bayes.multi_model_error import MultiModelError


class MultiAsymptoticModelError(MultiModelError):
    """
    Define a class that combines multiple asymptotic model errors
    """
    def __init__(self, yaml_file_metadata, yaml_file_data):
        """Create class combining different asymptotic model errors

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
        x_Mat_sensors = np.asarray(meta_data['x_Mat_sensors'])
#        mA = meta_data['mA']
#        mB = meta_data['mB']
        mA = np.asarray(meta_data['mA'])
        mB = np.asarray(meta_data['mB'])
        a = np.asarray(meta_data['a'])
        b = np.asarray(meta_data['b'])
        c = np.asarray(meta_data['c'])
        d = np.asarray(meta_data['d'])

        with open(yaml_file_data, "r") as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        data_MatB = np.asarray(data['MatB'])
        data_MatC = np.asarray(data['MatC'])
        data_Mat_B_C = np.asarray(data['Mat_B_C'])

        # assert(a == data['a'])
        # assert(b == data['b'])
        # assert(c == data['c'])
        # assert(d == data['d'])
        # assert(e == data['e'])

        single_model_error = AsymptoticModelError(x_Mat_sensors, data_MatB, data_MatC, data_Mat_B_C, a, b, c, d, mA, mB)  #, mA, mB
            # read in the parameters that are given in the experiment (in this case constant offset a)
        parameter = single_model_error.get_parameter_dict()
            # add parameters that are model parameters (not given in the experimental data)
        parameter.define("mA")
        parameter.define("mB") #! was b
        self.add(single_model_error, parameter)

        # define one shared latent parameter (only one parameter for all models)
        self.latent.add_by_name('mA')
        self.latent.add_by_name('mB') # was off. only one parameter before.


class TestOptimize(unittest.TestCase):
    def test_asymptotic_virtual_data(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_asymptotic_model_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiAsymptoticModelError('virtual_experiment_asymptotic_model_meta.yaml',
                                  'virtual_experiment_asymptotic_model_data.yaml')

        start_vector = np.array([0.98,-279.]) 
#, 0.98, -279 # one start value for mA in parameters check, list comprehension extended by steg for second parameter
# start vector must have 1 dimension
        result = least_squares(all_experiments_model_error, start_vector)

        # assert quadratic component of virtual experiment to be zero
        #self.assertAlmostEqual(d['mA'], 0.) # was there before.

        # check to obtain exact model coefficient
        self.assertAlmostEqual(result.x[0], d['mA'], 2) # was mA last value is accuracy in decimal places, was 5
        self.assertAlmostEqual(result.x[1], d['mB'], 2) # was mA last value is accuracy in decimal places, was 5
        # check cost function to be zero
        self.assertAlmostEqual(result.cost, 0)

    def test_asymptotic_virtual_data_split(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_asymptotic_model_with_noise_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiAsymptoticModelError('virtual_experiment_asymptotic_model_with_noise_meta.yaml',
                                  'virtual_experiment_asymptotic_model_with_noise_data.yaml')

        def weighted_all_experiments_model_error(prm_vector):
            all_model_errors = all_experiments_model_error.evaluate(prm_vector)
            #print('DEBUG prmvec',prm_vector)
            obj = 0.
            for key, me in all_model_errors.items():
                #print('DEBUG key',key)
                #print('DEBUG me',me)
                # obj = obj + np.dot(me['MatB'], me['MatB']) / d['sigma_noise_Mat'] ** 2 + \
               #       np.dot(me['MatC'], me['MatC']) / d['sigma_noise_Mat'] ** 2 \
                #         +  np.dot(me['Mat_B_C'], me['Mat_B_C']) / d['sigma_noise_Mat'] ** 2
                obj = obj + np.dot(me['Mat_B_C'], me['Mat_B_C']) / d['sigma_noise_Mat'] ** 2 
                #print('DEBUG obj',obj)
            return np.sqrt(obj)

        weighted_start_vector = np.array([0.98,-279.]) #was 0.98, -279# list extended by steg for second parameter, mA, mB
        weighted_result = least_squares(weighted_all_experiments_model_error, weighted_start_vector)

        # assert quadratic component of virtual experiment to be zero
        # self.assertAlmostEqual(d['a'], 0.) # no parameter to be zero here
        # check to obtain exact asymptotic coefficient
        self.assertAlmostEqual(weighted_result.x[0], d['mA'], 1) # was 4 in accuracy
        self.assertAlmostEqual(weighted_result.x[1], d['mB'], 1) # was 4 in accuracy in mA, mA or mB or both?


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
