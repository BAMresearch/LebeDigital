import logging
import numpy as np
from pathlib import Path
import os
import unittest
#from bayes.vb import *
import matplotlib.pyplot as plt
from theano_operation import LogLikeWithGrad
import pymc3 as pm
import theano.tensor as tt
from collections import OrderedDict
import copy
import matplotlib.cm as cm
import scipy.stats as st
import collections.abc
import yaml

"""
Todo
* include two types of sensors (function and derivative)
* this requires changing the structure of model error (remove model error)
    * should not return a pure vector,but
    * return forward_model[configuration][sensor]
    * return database[configuration][sensor]
    * assemble model error yourself
* include correlation between sensors 
* include correlation between configurations
* make sure to allow different configurations to have different sensors    
"""

class VirtualExperiment():
    def __init__(self, b, c, noise_std_function, correlation_length_function, noise_std_derivative,
                 correlation_length_derivative):
        self.b = b
        self.c = c
        self.noise_std_function = noise_std_function
        self.correlation_length_function = correlation_length_function
        self.noise_std_derivative = noise_std_derivative
        self.correlation_length_derivative = correlation_length_derivative

    def result(self, sensor_type, configuration, x):
        if sensor_type =="function":
            return self.b * x + configuration.parameters['a'] + self.c * x * x
        if sensor_type =="derivative":
            return 0 * x + self.b
        print("sensor type not known")
        exit(-1)

    def correlation_length(self, sensor_type):
        if sensor_type =="function":
            return self.correlation_length_function
        if sensor_type =="derivative":
            return self.correlation_length_derivative
        print("sensor type not known")
        exit(-1)

    def noise_std(self, sensor_type):
        if sensor_type == "function":
            return self.noise_std_function
        if sensor_type == "derivative":
            return self.noise_std_derivative
        print("sensor type not known")
        exit(-1)


class SensorGroup:
    def __init__(self, type, locations, correlation_length=None, tol=1e-4, ):
        self.type = type
        self.locations = locations
        if (len(locations)==0):
            print("zero number of sensors in group")
            exit(-1)
        self.correlation_length = correlation_length
        if correlation_length == None or correlation_length==0.:
            self.correlation = False
        else:
            self.correlation = True
            self.transformation = ModelErrorTransformation(locations, correlation_length, tol)

    def get_type(self):
        return self.type

    def is_correlated(self):
        return self.correlation


class Configuration:
    def __init__(self, a, sensor_groups):
        """
        For each individual a separate evaluation of the forward model is to be performed
        This could represent different time steps/transient computations,
        totally independent geometries, different (deterministic) boundary conditions
        Note that "a" only includes deterministic parameters
        The sensors are stored here as well, since different configurations might have
        totally different sensors
        """
        self.parameters = {} # deterministic parameters of the configuration, e.g. the time or given load
        self.parameters['a'] = a
        self.sensor_groups = sensor_groups

def correlation(locations, correlation_length):
    num_locations = locations.shape[0]
    c0 = np.repeat([locations], num_locations, axis=0)
    c1 = np.transpose(c0)
    r = c0 - c1
    # print("r",r)
    return np.exp(-r * r / (2. * correlation_length * correlation_length))

class ModelErrorTransformation:
    def __init__(self, locations, correlation_length, tol=1e-8):
        self.correlation = correlation(locations, correlation_length)
        self.tol = tol

        """
        Decompose the covariance matrix into its principal components
        Only keep the eigenvalues e with e > tol * largest eigenvalue

        Return the diagonal entries (representing the squares of the std_dev
        of independent random variables) and the corresponding eigenvectors  

        The full (correlated) sample vector X is then given by
        X = sum_i Phi_i * X_red,i with X_red,i being normal distributed with 
        zero mean and sigma2 given by the eigenvalues of the covariance matrix and Phi_i
        the corresponding eigenvalues
        """
        eigenvalues, eigen_vectors = np.linalg.eigh(self.correlation)
        treshold = tol * eigenvalues[-1]
        #("full eigenvalues", eigenvalues)
        # print("full eigenvectors", eigen_vectors)
        # print("threshold", treshold)
        self.reduced_eigenvalues = eigenvalues[eigenvalues > treshold]
        self.reduced_eigenvectors = eigen_vectors[:, eigenvalues > treshold]
        # print("eigenvectors", self.reduced_eigenvectors)
        print("reduced eigenvalues", self.reduced_eigenvalues)
        # scale the eigenvectors such that all reduced rv are unit variance
        # self.transformation = self.reduced_eigenvectors
        self.transformation = np.divide(self.reduced_eigenvectors, np.sqrt(self.reduced_eigenvalues))
        #print("Transformation", self.transformation)

    def dim_transformed_model_error(self):
        return len(self.reduced_eigenvalues)

    def num_noise_groups(self):
        return 1

    def original_cov(self, noise_post):
        '''
        computes the full covariance matrix of the sensor data
        COV = phiT * diag(variance of each eigenvalue) phi
        '''
        noise_variance_per_eigenvalue = self.noise_variance_per_eigenvalue(noise_post)
        return self.reduced_eigenvectors @ np.diag(noise_variance_per_eigenvalue) @ self.reduced_eigenvectors.T

    def original_var(self, noise_post):
        '''
        computes ony the diagonal entry of the full covariance matrix in self.original_cov
        '''
        noise_variance_per_eigenvalue = self.noise_variance_per_eigenvalue(noise_post)
        return np.sum(np.multiply(noise_variance_per_eigenvalue, np.square(self.reduced_eigenvectors)), axis=1)

    def noise_variance_per_eigenvalue(self, noise_post):
        '''
        returns (co)variance of the transformed model error (without scaling with eigenvalues, i.e. phiT*model_error)

        '''
        if len(noise_post.mean) > 1:
            print("noise posterior should only have a single entry")
            exit(-1)
        return 1. / noise_post.mean[0] * self.reduced_eigenvalues

    def transform_error(self, error):
        #print("error", error)
        #print("transformed error", error @ self.transformation)
        return error @ self.transformation

    def transform_jacobian(self, jacobian):
        return self.transformation.T @ jacobian


class Database:
    """
    is the interface to the database that stores the experimental data
    """
    def __init__(self, virtual_experiment):
        """
        this is now just adding virtual data, in reality there would be a database query
        and reading in the appropriate sensor information
        all the relevant meta data is stored either in the sensor object, e.g. the location of the
        sensor, or in the configuration  (global time steps, global load position, 
        or in our case the deterministic constant offset that varies between different configurations) 
        """
        self.virtual_experiment = virtual_experiment
        self.data= {}

    """ 
    returns a response for each configuration and sensor, as a 2D dict with the final result being an array 
    """
    def get_data(self, configurations):
        state = np.random.get_state()
        np.random.seed(42)
        #add noise (do not change the global random number generator)
        for configuration in configurations:
            if not configuration in self.data:
                self.data[configuration]= {}
                for sensor_group in configuration.sensor_groups:
                    self.data[configuration][sensor_group] = \
                        self.virtual_experiment.result(sensor_group.get_type(),configuration, sensor_group.locations)
                    correlation_length = self.virtual_experiment.correlation_length(sensor_group.get_type())
                    noise_std = self.virtual_experiment.noise_std(sensor_group.get_type())

                    if correlation_length==None or correlation_length<=0:
                        #uncorrelated noise
                        self.data[configuration][sensor_group] = self.data[configuration][sensor_group] \
                            +np.random.normal(0, noise_std, len(sensor_group.locations))
                    else:
                        fine_grid_locations = np.linspace(np.min(sensor_group.locations),
                                                          np.max(sensor_group.locations),
                                                          max(2,int((np.max(sensor_group.locations)-np.min(
                                                                sensor_group.locations))/correlation_length*50)))
                        noise_fine_grid = np.random.multivariate_normal(
                                np.zeros(len(fine_grid_locations)),
                                correlation(fine_grid_locations, correlation_length)*
                                    noise_std**2,  1)[0]
                        self.data[configuration][sensor_group] = self.data[configuration][sensor_group] + \
                            np.interp(sensor_group.locations,
                                           fine_grid_locations,
                                           noise_fine_grid)

        np.random.set_state(state)
        return self.data

class NumericalModel:
    """
    Calculates a model response for each sensor
    (configuration is fixed and given by the parameters in all_prms)
    """
    def __init__(self):
        pass

    def __call__(self, all_prms, sensor_groups):
        a = all_prms['a']
        b = all_prms['b']

        #Perform complex computation (not necessary here)

        #Helper function to extract sensor data
        def extract_sensor_value(sensor_type, x):
            if (sensor_type=="function"):
                return b * x + a
            if (sensor_type=="derivative"):
                return 0 * x + b
            print("Error in Helper function to extract sensor data")
            exit(-1)

        return {sensor_group : extract_sensor_value(sensor_group.get_type(), sensor_group.locations)  for
                sensor_group in
                sensor_groups}

class ModelError:
    def __init__(self, forward_model, database, configurations, configuration_correlation=None):
        """
        This calculates a model error for each configuration.
        """
        self.forward_model = forward_model
        self.database = database
        self.configurations = configurations
        # the data should be cached since it is not changing within the inference process
        self.data_of_all_sensors = self.database.get_data(self.configurations)
        self.configuration_correlation = configuration_correlation

    def __call__(self,latent_variables):
        model_error = {}
        for configuration in self.configurations:
            forward_me = self.forward_model(self._all_prms(configuration, latent_variables),configuration.sensor_groups)
            model_error[configuration] = {}
            for sensor_group in configuration.sensor_groups:
                model_error[configuration][sensor_group] = forward_me[sensor_group] - self.data_of_all_sensors[
                    configuration][sensor_group]
        return model_error

    def _all_prms(self, configuration, latent_prms):
        # In this example we assume that latent_parameters includes "b"
        # the configuration is providing "a"
        #
        # BUT: With the concept of _named_ parameters, this method must
        # include a change of data types
        #   all_prms:    A dictionary-like collection for convenient and
        #                more verbose access like all_prms.a
        #   latent_prms: A simple vector that is provided by pymc3 and VB
        #
        # This can be realized with a "ParameterView" and knowledge of which
        # of all_prms are latent - so basically the prior distribution.
        result = copy.deepcopy(configuration.parameters)
        result['b'] = latent_prms[0]
        return result

class VBModelError:
    def __init__(self, sensor_model_error):
        self.sensor_model_error = sensor_model_error

    def __call__(self, latent_prms):
        response = self.sensor_model_error(latent_prms)
        # just build a big vector with all the values
        values = np.array(list(response.values()))
        return values


class LogLike:
    """
    To be used within pymc3
    """

    def __init__(self, sensor_model_error):
        self.sme = sensor_model_error

    def __call__(self, prm):
        """
        VB requires a distribution of the  "precision" of the noise. To match
        that in the loglike, I modified equation
        
        return -0.5 * len(values) * np.log(2.0 * np.pi * (sigma ** 2)) - (
            0.5 * sigma ** 2
        ) * np.sum(np.square(values))

        with precision = 1. / sigma **2
        """
        noise_std_function = prm[-2]
        noise_std_derivative = prm[-1]
        latent_prms = prm[:-2]

        #[configuration][sensor]
        all_model_errors = self.sme(latent_prms)
        logLike = 0.
        if (self.sme.configuration_correlation==None):
            for configuration in self.sme.configurations:
                for sensor_group in configuration.sensor_groups:
                    if sensor_group.get_type() == "function":
                        noise_std = noise_std_function
                    else:
                        if sensor_group.get_type() == "derivative":
                            noise_std = noise_std_derivative
                        else:
                            print("Error in LogLike to extract sensor type")
                            exit(-1)
                    if (sensor_group.is_correlated()):
                        transformed_error = sensor_group.transformation.transform_error(
                            all_model_errors[configuration][sensor_group])
                        #print("transformed=",transformed_error)
                        logLike = logLike - 0.5 * len(transformed_error) * np.log(
                            2.0 * np.pi * (noise_std ** 2)) - \
                            (0.5 / noise_std ** 2) * np.sum(np.square(transformed_error))
                    else:
                        logLike = logLike - 0.5 * len(all_model_errors[configuration][sensor_group]) * np.log(2.0 * np.pi * (noise_std ** 2)) - \
                                  (0.5 / noise_std ** 2) * np.sum(np.square(all_model_errors[configuration][sensor_group]))
        else:
            print("correlation between configurations not yet implemented")
            exit(-1)
        return logLike


def run_pymc3(model_error):
    np.random.seed(182152)

    my_loglike = LogLike(model_error)
    logl = LogLikeWithGrad(my_loglike)
    with pm.Model():
        # Define priors !Make sure that this list corresponds to the right extraction for latent_parameters
        b = pm.Normal("b", 3.0, sd=1.0)
        noise_function_prior_mean = 0.2 #mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
        # good choice, then b=mean*(alpha-1), with alpha=4 this results in b=mean*3
        sigma_std_function = pm.InverseGamma("sigma_std_function", alpha=4., beta=noise_function_prior_mean*3)
        noise_derivative_prior_mean = 0.2 #mean=b/(a-1) var=b**2/(a-1)**/a -> a=mean**2/var but>2, so usually a=4 is a
        sigma_std_derivative = pm.InverseGamma("sigma_std_derivative", alpha=4., beta=noise_derivative_prior_mean*3)

        theta = tt.as_tensor_variable([b, sigma_std_function, sigma_std_derivative])
#        pm.DensityDist("likelihood", lambda v: logl(v), observed={"v": theta})
        pm.Potential("likelihood", logl(theta))

        # Inference!
        trace = pm.sample(
            draws=2000,
            step=pm.Metropolis(),
            chains=4,
            tune=100,
            discard_tuned_samples=True,
        )

    print(trace)
    s = pm.summary(trace)
    print(s)
    means = s["mean"]
    sds = s["sd"]
    print(f"PM: Posterior for 'b' = {means['b']:6.3f} with sd {sds['b']:6.3f}.")
    print(f"PM: Posterior for 'sigma_std_function' = {means['sigma_std_function']:6.3f} with sd "
          f"{sds['sigma_std_function']:6.3f}.")
    print(f"PM: Posterior for 'sigma_std_derivative' = {means['sigma_std_derivative']:6.3f} with sd "
          f"{sds['sigma_std_derivative']:6.3f}.")
    #print(trace.stat_names)
    #accept = trace.get_sampler_stats('accept')
    #print("accept", accept)
    #pm.traceplot(trace, priors=[b.distribution,sigma_std_function.distribution, sigma_std_derivative.distribution]);
    plt.show()
    return trace

def plot_inv_gamma():
    print("mean", st.invgamma(2, scale=1e-9).mean())
    x = np.linspace(0.0, 1, 1000)
    fig, ax = plt.subplots()
    f = lambda alpha, beta: st.invgamma.pdf(x, alpha, scale=beta)
    plot_pdf = lambda alpha, beta: ax.plot(x, f(alpha, beta),
                                           label=r'$\alpha$={0}, $\beta$={1}'.format(alpha, beta))
    plot_pdf(2, 0.1)
    # plot_pdf(1., 1e-2)
    # plot_pdf(1.1, 1e-2)
    plt.legend(loc='upper right', frameon=False)
    ax.set(xlim=[0, 1], xlabel='x', ylabel='f(x)')
    plt.show()

def setup(options):
    # Define some sensors
    x_function_sensors = np.linspace(0, 1, options['n_function_sensors'])
    x_derivative_sensors = np.linspace(0, 1, options['n_derivative_sensors'])

    sensor_groups =[]
    if options['n_function_sensors']>0:
        sensor_groups.append(SensorGroup("function"  ,x_function_sensors   , correlation_length=options['correlation_length_function']))
    if options['n_derivative_sensors']>0:
        sensor_groups.append(SensorGroup("derivative", x_derivative_sensors, correlation_length=options['correlation_length_derivative']))

    # Our model can evaluate all sensors
    forward_model = NumericalModel()

    #in this implementation there is no database, but the database functions
    # just compute the virtual data based on a virtual experiment
    virtual_experiment = VirtualExperiment(
        b=options['virtual_experiment_b'],
        c=options['virtual_experiment_c'],
        noise_std_function=options['virtual_experiment_noise_std_function'],
        correlation_length_function=options['virtual_experiment_correlation_length_function'],
        noise_std_derivative = options['virtual_experiment_noise_std_derivative'],
        correlation_length_derivative=options['virtual_experiment_correlation_length_derivative'])
    data_base = Database(virtual_experiment)

    # define configurations (in our case offset a=0 and a=10, could be different load configurations
    # or times
    configurations = [Configuration(parameter_a,sensor_groups) for parameter_a in np.linspace(0,4,options['n_configurations'])]

    model_error = ModelError(forward_model, data_base, configurations)

    return [model_error, virtual_experiment, configurations, forward_model, data_base]

class Opts(collections.Mapping):
    """
    Dictionary-like gimmick class that, in contrast to a regular dictionary,:
        * is always connected to a parameter yaml file
        * allows you to access the content of the yaml parameter file
          via:
              * dict-like access      -- options["my_prm"]
              * attribute-like access -- options.my_prm
        * allows you to access the parameters case-insensitive
        * explicitly forbids modifying its parameters
    """

    def __init__(self, prm_file):
        import yaml

        self._data = {}
        with open(prm_file, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)

        for key, value in d.items():
            self._data[key.lower()] = value

    def __getitem__(self, key):
        return self._data[key.lower()]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __getattr__(self, name):
        return self[name]

    def __setitem__(*args):
        raise RuntimeError("You are not allowed to modify options.")

    def value_or(self, key, default):
        try:
            return self._data[key]
        except:
            return default