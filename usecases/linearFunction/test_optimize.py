import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import yaml
from linear_model_error import LinearModelError
from bayes.multi_model_error import MultiModelError
from theano_operation import LogLikeWithGrad
import pymc3 as pm
import theano.tensor as tt


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
    def test_linear(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_linear_meta.yaml',
                                  'virtual_experiment_linear_data.yaml')

        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        # assert quadratic component of virtual experiment to be zero
        self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        self.assertAlmostEqual(result.x[0], d['b'], 7)
        # check cost function to be zero
        self.assertAlmostEqual(result.cost, 0)

    def test_linear_with_noise_optimize(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_with_noise_optimize_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_linear_with_noise_optimize_meta.yaml',
                                  'virtual_experiment_linear_with_noise_optimize_data.yaml')

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
        self.assertAlmostEqual(weighted_result.x[0], d['b'], 4)

    def test_linear_with_noise_mcmc(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_with_noise_mcmc_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_linear_with_noise_mcmc_meta.yaml',
                                  'virtual_experiment_linear_with_noise_mcmc_data.yaml')

        # add noise term as local parameter to all model errors
        #todo change latent in multimodelerror, since this is how we import the complete list
        #todo do we have to add hyperparameters to the model parameters, since they are not part of the forward model
        #howto access the parameter list outside (e.g. to add a model parameter)
        for key, param_list in all_experiments_model_error.latent._all_parameters.items():
            param_list.define('sigma_noise_f', 1.)
            param_list.define('sigma_noise_df', 1.)

        # add noise term as global shared parameter (all model errors have the same)
        all_experiments_model_error.latent.add_by_name('sigma_noise_f')
        all_experiments_model_error.latent.add_by_name('sigma_noise_df')

        class LogLike:
            """
            To be used within pymc3
            """
            def __init__(self, multi_model_error):
                self.multi_model_error = multi_model_error

            def __call__(self, prm_vector):
                """
                return -0.5 * len(values) * np.log(2.0 * np.pi * (sigma ** 2)) - (
                    0.5 * sigma ** 2
                ) * np.sum(np.square(values))

                with precision = 1. / sigma **2
                """
                # should evaluate return a copy or update the parameter vector ? Standard is now a copy
                # but this is not returned and thus cannot be used in the evaluation of loglike
                multi_model_error_result = self.multi_model_error.evaluate(prm_vector)
                log_like = 0.
                for key, me in multi_model_error_result.items():
                    param = all_experiments_model_error.latent._all_parameters[key]

                    log_like = log_like - 0.5 * (
                            len(me['f']) * np.log(2.0 * np.pi * (param['sigma_noise_f'] ** 2)) \
                            + np.sum(np.square(me['f']/param['sigma_noise_f'] ** 2))
                            + len(me['df']) * np.log(2.0 * np.pi * (param['sigma_noise_df'] ** 2)) \
                            + np.sum(np.square(me['df']/param['sigma_noise_df'] ** 2))
                            )
                return log_like

        np.random.seed(42)

        my_loglike = LogLike(all_experiments_model_error)
        logl = LogLikeWithGrad(my_loglike)
        with pm.Model():
            # Define priors !Make sure that this list corresponds to the right extraction for latent_parameters
            b = pm.Normal("b", 5.0, sd=1.0)
            noise_function_prior_mean = 5.  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
            # good choice, then b=mean*(alpha-1), with alpha=4 this results in b=mean*3
            sigma_noise_function = pm.InverseGamma("sigma_noise_function", alpha=4., beta=noise_function_prior_mean * 3)
            noise_derivative_prior_mean = 0.5  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
            sigma_noise_derivative = pm.InverseGamma("sigma_noise_derivative", alpha=4.,
                                                   beta=noise_derivative_prior_mean * 3)

            theta = tt.as_tensor_variable([b, sigma_noise_function, sigma_noise_derivative])
            pm.Potential("likelihood", logl(theta))

            # Inference!
            trace = pm.sample(
                draws=2000,
                step=pm.Metropolis(),
                chains=4,
                tune=500, #used for tuning the proposal density to get specified acceptance rate (65?)
                discard_tuned_samples=True,
            )

        print(trace)
        s = pm.summary(trace)
        print(s)
        means = s["mean"]
        sds = s["sd"]
        print(f"PM: Posterior for 'b' = {means['b']:6.3f} with sd {sds['b']:6.3f}.")
        print(f"PM: Posterior for 'sigma_std_function' = {means['sigma_noise_function']:6.3f} with sd "
              f"{sds['sigma_noise_function']:6.3f}.")
        print(f"PM: Posterior for 'sigma_std_derivative' = {means['sigma_noise_derivative']:6.3f} with sd "
              f"{sds['sigma_noise_derivative']:6.3f}.")
        # print(trace.stat_names)
        # accept = trace.get_sampler_stats('accept')
        # print("accept", accept)
        #pm.traceplot(trace, priors=[b.distribution,sigma_noise_function.distribution, sigma_noise_derivative.distribution]);
        #plt.show()
        #return trace

        # assert quadratic component of virtual experiment to be zero
        # self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        # self.assertAlmostEqual(weighted_result.x[0], d['b'], 4)

    def test_quadratic(self):
        all_experiments_model_error = \
            MultiLinearModelError('virtual_experiment_quadratic_meta.yaml',
                                  'virtual_experiment_quadratic_data.yaml')
        start_vector = np.array([0.7])
        result = least_squares(all_experiments_model_error, start_vector)

        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_quadratic_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        # check to obtain linear coefficient using sympy
        # solve(diff(integrate(((c*x**2+b*x)-bbar*x)**2,(x,0,1)),bbar),bbar)
        # accuracy is limited due to integration error with a limited number of samples
        self.assertAlmostEqual(result.x[0], d['b']+3./4.*d['c'], 6)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    unittest.main()
