import unittest
import numpy as np
import yaml

class VirtualLinearModelExperiment:
    """
    create virtual test data for a function a+b*x+c*x**2
    evalute the function at positions x_function
    evaluate the derivate at positions x_derivatives
    a is a variable that is different in each experiment and
    creates a constant offset (different, but known)
    """
    def __init__(self, b, c, num_function_sensors, num_derivative_sensors, center=False):
        self.b = b
        self.c = c
        if center is False:
            self.x_function = np.linspace(0, 1, num_function_sensors)
            self.x_derivative = np.linspace(0, 1, num_derivative_sensors)
        else:
            if num_function_sensors>0:
                h_function = 1. / num_function_sensors
                self.x_function = np.linspace(0, 1, num_function_sensors, endpoint=False) + 0.5 * h_function
            else:
                self.x_function = np.empty( shape=(0, 0))

            if num_derivative_sensors>0:
                h_derivative = 1. / num_derivative_sensors
                self.x_derivative = np.linspace(0, 1, num_derivative_sensors, endpoint=False) + 0.5 * h_derivative
            else:
                self.x_derivative = np.empty( shape=(0, 0))

    def write_data_to_yaml(self, virtual_experimental_data_file, a):
        f_x = a * np.ones(len(self.x_function)) + \
               self.b * self.x_function + self.c * np.square(self.x_function)
        df_x = self.b * np.ones(len(self.x_derivative)) + 2. * self.c * self.x_derivative
        data = {
            "x_function": self.x_function.tolist(),
            "x_derivative": self.x_derivative.tolist(),
            "a": a,
            "f": f_x.tolist(),
            "df": df_x.tolist()
        }
        with open(virtual_experimental_data_file, "w") as f:
            d = yaml.dump(data, f, default_flow_style=None)
        f.close()

    def write_model_to_yaml(self, virtual_experiment_model_file):
        data = {
            "b": self.b,
            "c": self.c
        }
        with open(virtual_experiment_model_file, "w") as f:
            d = yaml.dump(data, f, default_flow_style=None)
        f.close()

class TestVirtualExperiments(unittest.TestCase):

    def test_linear_virtual_data(self):
        virtual_experiment = VirtualLinearModelExperiment(b=3, c=0,
            num_function_sensors=10, num_derivative_sensors=6)
        virtual_experiment.write_model_to_yaml(f"virtual_linear_experiment_model.yaml")

        #generate 10 experiments with different offset a (here assumed to be 2*number)
        for experiment in range(10):
            virtual_experiment.write_data_to_yaml(f"virtual_linear_experiment_data_{experiment}.yaml", a=2*experiment)

    def test_quadratic_virtual_data(self):
        virtual_experiment = VirtualLinearModelExperiment(b=0, c=2, num_function_sensors=1000,
                                                          num_derivative_sensors=0, center=True)
        virtual_experiment.write_model_to_yaml(f"virtual_quadratic_experiment_model.yaml")

        # generate 10 experiments with different offset a (here assumed to be 2*number)

        for experiment in range(1):
            virtual_experiment.write_data_to_yaml(f"virtual_quadratic_experiment_data_{experiment}.yaml", a=0 * experiment)


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)
    unittest.main()
