import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import yaml
from linear_model_error import LinearModelError
from bayes.inference_problem import InferenceProblem
from theano_operation import LogLikeWithGrad
import pymc3 as pm
import theano.tensor as tt


def create_inference_problem(yaml_file_metadata, yaml_file_data):
    with open(yaml_file_metadata, "r") as f:
        meta_data = yaml.load(f, Loader=yaml.FullLoader)

    x_function = np.asarray(meta_data['x_function'])
    x_derivative = np.asarray(meta_data['x_derivative'])
    all_a = np.asarray(meta_data['all_a'])

    with open(yaml_file_data, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    problem = InferenceProblem()
    for idx, a in enumerate(all_a):
        data_f = np.asarray(data[idx]['f'])
        data_df = np.asarray(data[idx]['df'])
        assert(a == data[idx]['a'])
        single_model_error = LinearModelError(x_function, x_derivative, data_f, data_df, a)
        problem.add_model_error(single_model_error, key=idx)
        # set the parameters that are given in the experiment (in this case constant offset a)
        parameter_list = single_model_error.parameter_list
        parameter_list['a'] = a
        # add the latent parameters to the joint latent list
        #problem.latent['b'].add(single_model_error, 'b')
    #set b as a global latent variable for all model errors that have a parameter b
    problem.define_shared_latent_parameter_by_name("b")
    return problem

def latent_as_tensor_variable(self, model, prior_pymc3):
    with model:
        #array = [prior_pymc3[key] for key, var in prior_pymc3.items()]
        array = []
        for latent_parameters in self.latent._index_mapping:
            # latent_parameters is now a dict of all parameters with key related to the model key and value being the name of the variable in this list
            # get value of first iterator, which is currently supposed to be all identical
            # TODO we need to define for each latent variable a new name that is then used to identify the prior, or set the initial value
            var_name = next(iter(latent_parameters.values()))
            array.append(prior_pymc3[var_name])
        return tt.as_tensor_variable(array)

class ConcatenatedME:
    def __init__(self, inference_problem):
        self.inference_problem = inference_problem

    def __call__(self, parameter_vector):
        me = self.inference_problem(parameter_vector)
        tmp = np.concatenate([sensor_me for experiment in me.values() for sensor_me in experiment.values()])
        return tmp

class TestOptimize(unittest.TestCase):
    def test_linear(self):
        virtual_experiment_metadata_file = \
            Path(__file__).parent / 'virtual_experiment_linear_meta.yaml'
        with open(virtual_experiment_metadata_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        inference_problem = create_inference_problem('virtual_experiment_linear_meta.yaml',
                                  'virtual_experiment_linear_data.yaml')

        # #set starting values, maybe move this to the library inference_problem.latent['b'].set(0.7)
        # latent = inference_problem.latent['b']
        # for parameter_list, name in latent:
        #     parameter_list[name] = 0.7
        # #extract the current values

        #extract the length of the complete latent parameter vector by looping over all latenten parameters and
        # summing up their length (a parameter might be a vector)
        latent_length = sum(latent_par.N for latent_par in inference_problem.latent.values())
        #init latent vector (not zero, so when forgetting to initialize a value, an error is returned
        start_vector = np.empty(latent_length)
        #fill the the global vector with the inital value of the variable (b) based on its global position
        start_vector[inference_problem.latent['b'].global_index_range()] = 0.7

        concatenated_me = ConcatenatedME(inference_problem)
        result = least_squares(concatenated_me, start_vector)

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

        inference_problem = create_inference_problem('virtual_experiment_linear_with_noise_optimize_meta.yaml',
                                  'virtual_experiment_linear_with_noise_optimize_data.yaml')

        def weighted_all_experiments_model_error(prm_vector):
            all_model_errors = inference_problem(prm_vector)
            obj = 0.
            for me in all_model_errors.values():
                obj = obj + np.dot(me['f'], me['f']) / d['sigma_noise_function'] ** 2 + \
                      np.dot(me['df'], me['df']) / d['sigma_noise_derivative'] ** 2
            return np.sqrt(obj)

        #extract the length of the complete latent parameter vector by looping over all latenten parameters and
        # summing up their length (a parameter might be a vector)
        latent_length = sum(latent_par.N for latent_par in inference_problem.latent.values())
        #init latent vector (not zero, so when forgetting to initialize a value, an error is returned
        start_vector = np.empty(latent_length)
        #fill the the global vector with the inital value of the variable (b) based on its global position
        start_vector[inference_problem.latent['b'].global_index_range()] = 0.7

        #optimize
        weighted_result = least_squares(weighted_all_experiments_model_error, start_vector)

        # assert quadratic component of virtual experiment to be zero
        self.assertAlmostEqual(d['c'], 0.)
        # check to obtain exact linear coefficient
        self.assertAlmostEqual(weighted_result.x[0], d['b'], 4)

   #  def test_linear_with_noise_mcmc(self):
   #      virtual_experiment_metadata_file = \
   #          Path(__file__).parent / 'virtual_experiment_linear_with_noise_mcmc_meta.yaml'
   #      with open(virtual_experiment_metadata_file, "r") as f:
   #          d = yaml.load(f, Loader=yaml.FullLoader)
   #
   #      inference_problem = create_inference_problem('virtual_experiment_linear_with_noise_mcmc_meta.yaml',
   #                                'virtual_experiment_linear_with_noise_mcmc_data.yaml')
   #
   #      # add noise term as local parameter to all model errors
   #      # currently we have to add hyperparameters to the model parameters, since they are not part of the forward model
   #      for param_list in all_experiments_model_error.linearModelParameters:
   #          param_list.define('sigma_noise_f', 1.)
   #          param_list.define('sigma_noise_df', 1.)
   #
   #      #define log likelihood
   #      class LogLike:
   #          """
   #          To be used within pymc3
   #          """
   #          def __init__(self, multi_model_error):
   #              self.multi_model_error = multi_model_error
   #
   #          def __call__(self, prm_vector):
   #              """
   #              return -0.5 * len(values) * np.log(2.0 * np.pi * (sigma ** 2)) - (
   #                  0.5 * sigma ** 2
   #              ) * np.sum(np.square(values))
   #
   #              with precision = 1. / sigma **2
   #              """
   #
   #              multi_model_error_result = self.multi_model_error.evaluate(prm_vector)
   #              log_like = 0.
   #              for key, me in multi_model_error_result.items():
   #                  param = all_experiments_model_error.latent._all_parameters[key]
   #
   #                  log_like = log_like - 0.5 * (
   #                          len(me['f']) * np.log(2.0 * np.pi * (param['sigma_noise_f'] ** 2)) \
   #                          + np.sum(np.square(me['f']/param['sigma_noise_f'] ** 2))
   #                          + len(me['df']) * np.log(2.0 * np.pi * (param['sigma_noise_df'] ** 2)) \
   #                          + np.sum(np.square(me['df']/param['sigma_noise_df'] ** 2))
   #                          )
   #              return log_like
   #
   #      np.random.seed(42)
   #
   #      #instantiate log likelihood with the model error
   #      loglike_without_gradient = LogLike(all_experiments_model_error)
   #      #generate a likelihood function with an automatic generation of the gradient using finite differences
   #      loglike = LogLikeWithGrad(loglike_without_gradient)
   #      pyMC3_prior = [None] * 3 # inititalize array of pyMC3 priors with the number of latent parameters
   #      # add noise term as global shared parameter (all model errors have the same)
   #      with pm.Model() as model:
   #          # Define priors, add_by_name returns the index in the global parameter vector
   #          # this ensures that the pymc3 priors are ordered in the same way as the global vector
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('b')] = pm.Normal("b", 5.0, sd=1.0)
   #          noise_function_prior_mean = 5.  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
   #          # good choice, then b=mean*(alpha-1), with alpha=4 this results in b=mean*3
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('sigma_noise_f')] = pm.InverseGamma("sigma_noise_f", alpha=4., beta=noise_function_prior_mean * 3)
   #          noise_derivative_prior_mean = 0.5  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('sigma_noise_df')] = pm.InverseGamma("sigma_noise_df", alpha=4.,
   #                                                 beta=noise_derivative_prior_mean * 3)
   #          theta = tt.as_tensor_variable(pyMC3_prior)
   #          pm.Potential("likelihood", loglike(theta))
   #
   #          # Inference!
   #          trace = pm.sample(
   #              draws=1000,
   #              step=pm.Metropolis(),
   #              chains=4,
   #              tune=500, #used for tuning the proposal density to get specified acceptance rate (65%?)
   #              discard_tuned_samples=True,
   #          )
   #
   #      print(trace)
   #      s = pm.summary(trace)
   #      print(s)
   #      means = s["mean"]
   #      sds = s["sd"]
   #      print(f"PM: Posterior for 'b' = {means['b']:6.3f} with sd {sds['b']:6.3f}.")
   #      print(f"PM: Posterior for 'sigma_std_f' = {means['sigma_noise_f']:6.3f} with sd "
   #            f"{sds['sigma_noise_f']:6.3f}.")
   #      print(f"PM: Posterior for 'sigma_std_df' = {means['sigma_noise_df']:6.3f} with sd "
   #            f"{sds['sigma_noise_df']:6.3f}.")
   #      # print(trace.stat_names)
   #      # accept = trace.get_sampler_stats('accept')
   #      # print("accept", accept)
   #      #pm.traceplot(trace, priors=[b.distribution,sigma_noise_function.distribution, sigma_noise_derivative.distribution]);
   #      #plt.show()
   #      #return trace
   #
   #      # assert quadratic component of virtual experiment to be zero
   #      # self.assertAlmostEqual(d['c'], 0.)
   #      # check to obtain exact linear coefficient
   #      # self.assertAlmostEqual(weighted_result.x[0], d['b'], 4)
   # # def test_linear_with_noise_mcmc(self):
   #      virtual_experiment_metadata_file = \
   #          Path(__file__).parent / 'virtual_experiment_linear_with_noise_mcmc_meta.yaml'
   #      with open(virtual_experiment_metadata_file, "r") as f:
   #          d = yaml.load(f, Loader=yaml.FullLoader)
   #
   #      inference_problem = create_inference_problem('virtual_experiment_linear_with_noise_mcmc_meta.yaml',
   #                                'virtual_experiment_linear_with_noise_mcmc_data.yaml')
   #
   #      # add noise term as local parameter to all model errors
   #      # currently we have to add hyperparameters to the model parameters, since they are not part of the forward model
   #      for param_list in all_experiments_model_error.linearModelParameters:
   #          param_list.define('sigma_noise_f', 1.)
   #          param_list.define('sigma_noise_df', 1.)
   #
   #      #define log likelihood
   #      class LogLike:
   #          """
   #          To be used within pymc3
   #          """
   #          def __init__(self, multi_model_error):
   #              self.multi_model_error = multi_model_error
   #
   #          def __call__(self, prm_vector):
   #              """
   #              return -0.5 * len(values) * np.log(2.0 * np.pi * (sigma ** 2)) - (
   #                  0.5 * sigma ** 2
   #              ) * np.sum(np.square(values))
   #
   #              with precision = 1. / sigma **2
   #              """
   #
   #              multi_model_error_result = self.multi_model_error.evaluate(prm_vector)
   #              log_like = 0.
   #              for key, me in multi_model_error_result.items():
   #                  param = all_experiments_model_error.latent._all_parameters[key]
   #
   #                  log_like = log_like - 0.5 * (
   #                          len(me['f']) * np.log(2.0 * np.pi * (param['sigma_noise_f'] ** 2)) \
   #                          + np.sum(np.square(me['f']/param['sigma_noise_f'] ** 2))
   #                          + len(me['df']) * np.log(2.0 * np.pi * (param['sigma_noise_df'] ** 2)) \
   #                          + np.sum(np.square(me['df']/param['sigma_noise_df'] ** 2))
   #                          )
   #              return log_like
   #
   #      np.random.seed(42)
   #
   #      #instantiate log likelihood with the model error
   #      loglike_without_gradient = LogLike(all_experiments_model_error)
   #      #generate a likelihood function with an automatic generation of the gradient using finite differences
   #      loglike = LogLikeWithGrad(loglike_without_gradient)
   #      pyMC3_prior = [None] * 3 # inititalize array of pyMC3 priors with the number of latent parameters
   #      # add noise term as global shared parameter (all model errors have the same)
   #      with pm.Model() as model:
   #          # Define priors, add_by_name returns the index in the global parameter vector
   #          # this ensures that the pymc3 priors are ordered in the same way as the global vector
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('b')] = pm.Normal("b", 5.0, sd=1.0)
   #          noise_function_prior_mean = 5.  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
   #          # good choice, then b=mean*(alpha-1), with alpha=4 this results in b=mean*3
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('sigma_noise_f')] = pm.InverseGamma("sigma_noise_f", alpha=4., beta=noise_function_prior_mean * 3)
   #          noise_derivative_prior_mean = 0.5  # mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
   #          pyMC3_prior[all_experiments_model_error.latent.add_by_name('sigma_noise_df')] = pm.InverseGamma("sigma_noise_df", alpha=4.,
   #                                                 beta=noise_derivative_prior_mean * 3)
   #          theta = tt.as_tensor_variable(pyMC3_prior)
   #          pm.Potential("likelihood", loglike(theta))
   #
   #          # Inference!
   #          trace = pm.sample(
   #              draws=1000,
   #              step=pm.Metropolis(),
   #              chains=4,
   #              tune=500, #used for tuning the proposal density to get specified acceptance rate (65%?)
   #              discard_tuned_samples=True,
   #          )
   #
   #      print(trace)
   #      s = pm.summary(trace)
   #      print(s)
   #      means = s["mean"]
   #      sds = s["sd"]
   #      print(f"PM: Posterior for 'b' = {means['b']:6.3f} with sd {sds['b']:6.3f}.")
   #      print(f"PM: Posterior for 'sigma_std_f' = {means['sigma_noise_f']:6.3f} with sd "
   #            f"{sds['sigma_noise_f']:6.3f}.")
   #      print(f"PM: Posterior for 'sigma_std_df' = {means['sigma_noise_df']:6.3f} with sd "
   #            f"{sds['sigma_noise_df']:6.3f}.")
   #      # print(trace.stat_names)
   #      # accept = trace.get_sampler_stats('accept')
   #      # print("accept", accept)
   #      #pm.traceplot(trace, priors=[b.distribution,sigma_noise_function.distribution, sigma_noise_derivative.distribution]);
   #      #plt.show()
   #      #return trace
   #
   #      # assert quadratic component of virtual experiment to be zero
   #      # self.assertAlmostEqual(d['c'], 0.)
   #      # check to obtain exact linear coefficient
   #      # self.assertAlmostEqual(weighted_result.x[0], d['b'], 4)

    def test_quadratic(self):
        inference_problem = create_inference_problem('virtual_experiment_quadratic_meta.yaml',
                                  'virtual_experiment_quadratic_data.yaml')

        #extract the length of the complete latent parameter vector by looping over all latenten parameters and
        # summing up their length (a parameter might be a vector)
        latent_length = sum(latent_par.N for latent_par in inference_problem.latent.values())
        #init latent vector (not zero, so when forgetting to initialize a value, an error is returned
        start_vector = np.empty(latent_length)
        #fill the the global vector with the inital value of the variable (b) based on its global position
        start_vector[inference_problem.latent['b'].global_index_range()] = 0.7

        concatenated_me = ConcatenatedME(inference_problem)
        result = least_squares(concatenated_me, start_vector)

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
