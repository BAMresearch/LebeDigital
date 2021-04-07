import unittest
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy.optimize import newton
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
       # mA = np.asarray(meta_data['mA'])
       # mB = np.asarray(meta_data['mB'])
        a = np.asarray(meta_data['a'])
        b = np.asarray(meta_data['b'])
        c = np.asarray(meta_data['c'])
        d = np.asarray(meta_data['d'])
        e = np.asarray(meta_data['e'])

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
######################################################
        single_model_error = AsymptoticModelError(x_Mat_sensors, data_MatB, data_MatC, data_Mat_B_C, a, b, c, d, e) #, mA, mB)  #, mA, mB
######################################################
            # read in the parameters that are given in the experiment (in this case a, b, c, d, e)
        parameter = single_model_error.get_parameter_dict()
            # add parameters that are model parameters (not given in the experimental data)
        parameter.define("mA") # this correct
        parameter.define("mB") # this correct
        self.add(single_model_error, parameter)

        # define one shared latent parameter (only one parameter b for all linear models)
        self.latent.add_by_name('mA') # this is necessary
        self.latent.add_by_name('mB') # was off. only one parameter before.


class TestOptimize(unittest.TestCase):
    def test_asymptotic_virtual_data(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_asymptotic_model_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            dat = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = MultiAsymptoticModelError('virtual_experiment_asymptotic_model_meta.yaml', 'virtual_experiment_asymptotic_model_data.yaml')

        start_vector = np.array([ 0.5, -279.])    #, 0.98, -279 #  list comprehension extended by steg for second parameter
        result = least_squares(all_experiments_model_error, start_vector, method='trf',verbose=0) # method trf
        # method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08, x_scale='jac', diff_step=None,
        # result =  newton(all_experiments_model_error, start_vector)
        # result = curve_fit(all_experiments_model_error, start_vector)
        #print('DEBUG single model',result)

        ####### boundary conditions. accuracy criteria #######
        # assert quadratic component of virtual experiment to be zero
        # self.assertAlmostEqual(d['c'], 0.)

        # check to obtain exact model coefficient
        #self.assertAlmostEqual(result.x[0], dat['mA'], 4) # was mA last value is accuracy in decimal places, was 5
        #self.assertAlmostEqual(result.x[1], dat['mB'], 4) # was mA last value is accuracy in decimal places, was 5
        # check cost function to be zero
        self.assertAlmostEqual(result.cost, 0., 2)

    def test_asymptotic_virtual_data_split(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_asymptotic_model_with_noise_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            dat = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = MultiAsymptoticModelError('virtual_experiment_asymptotic_model_with_noise_meta.yaml', 'virtual_experiment_asymptotic_model_with_noise_data.yaml')

        def weighted_all_experiments_model_error(prm_vector):
            all_model_errors = all_experiments_model_error.evaluate(prm_vector)

            obj = 0.
            for key, me in all_model_errors.items():
                obj = obj + np.dot(me['Mat_B_C'], me['Mat_B_C']) / dat['sigma_noise_Mat'] ** 2
                #print('DEBUG obj',obj)
            return np.sqrt(obj)

        weighted_start_vector = np.array([0.98, -279.]) #was 0.98, -279# list extended by steg for second parameter, mA, mB
        weighted_result = least_squares(weighted_all_experiments_model_error, weighted_start_vector, method='trf', x_scale='jac', tr_solver='lsmr')
        # method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08, x_scale='jac', diff_step=None, max_nfev=10000

        # assert quadratic component of virtual experiment to be zero
        # self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact asymptotic coefficient
 #       self.assertAlmostEqual(weighted_result.x[0], dat['mA'], 4) # was 4 in accuracy
#        self.assertAlmostEqual(weighted_result.x[1], dat['mB'], 4) # was 4 in accuracy in mA, mA or mB or both?
        #self.assertAlmostEqual(weighted_result.cost, 0., 2)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
