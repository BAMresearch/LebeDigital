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
    def __init__(self, b, c, num_function_sensors, num_derivative_sensors):
        self.b = b
        self.c = c
        self.x_function = np.linspace(0, 1, num_function_sensors)
        self.x_derivative = np.linspace(0, 1, num_derivative_sensors)

    def write_to_yaml(self, virtual_experimental_data_file, a):
        f_x = a * np.ones(len(self.x_function)) + \
               self.b * self.x_function
        df_x = self.b * np.ones(len(self.x_derivative))
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


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)

    virtual_experiment = VirtualLinearModelExperiment(b=3, c=0,
        num_function_sensors=10, num_derivative_sensors=6)

    #generate 10 experiments with different offset a (here assumed to be 2*number)
    for experiment in range(10):
        virtual_experiment.write_to_yaml(f"LinearModelExperiment_{experiment}.yaml", a=2*experiment)

