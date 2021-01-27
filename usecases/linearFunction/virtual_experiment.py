import numpy as np
import yaml


class VirtualLinearModelExperiment:
    """create virtual test data for a function a+b*x+c*x**2

    Note:
        The model is actually nonlinear due to the quadratic component c that represents a model bias

    Attributes:
        all_a (np.array): array with offsets (for each entry in a, there is a corresponding measurement
            performed at x_function and x_derivative)
        b (float): linear coefficient of the model
        c (float): quadratic coefficient of the model
        x_function(np.array): positions of function sensors in the interval [0,1]
        x_derivative(np.array): positions of derivative sensors in the interval [0,1]
    """

    def __init__(self, virtual_experiment_metadata_yaml):
        """Create virtual experiment
        Args:
            virtual_experiment_metadata_yaml: meta data file to generate the data from
        """
        with open(virtual_experiment_metadata_yaml, "r") as f:
            d = yaml.load(f, Loader=yaml.FullLoader)
        self.all_a = np.asarray(d['all_a'])
        self.b = d['b']
        self.c = d['c']
        self.x_function = np.asarray(d['x_function'])
        self.x_derivative = np.asarray(d['x_derivative'])

    def write_data_to_yaml(self, virtual_experimental_data_file):
        """Write virtual sensor data to yaml file

        Args:
            virtual_experimental_data_file(str): file location

        Returns:
            writes the experimental data to the yaml file
        """
        with open(virtual_experimental_data_file, "w") as f:
            data = {}
            for index, a in enumerate(self.all_a):
                f_x = a * np.ones(len(self.x_function)) \
                    + self.b * self.x_function + self.c * np.square(self.x_function)
                df_x = self.b * np.ones(len(self.x_derivative)) + 2. * self.c * self.x_derivative
                data[index] = {
                    "a": float(a),
                    "f": f_x.tolist(),
                    "df": df_x.tolist()
                }
            yaml.dump(data, f, default_flow_style=None)


def main():
    # create experiment with an exactly linear model (no model bias)
    virtual_experiment = VirtualLinearModelExperiment("virtual_experiment_linear_model_meta.yaml")
    # write experiments raw data to file
    virtual_experiment.write_data_to_yaml("virtual_experiment_linear_model_data.yaml")

    # create data with an quadratic model (thus the linear model has a model bias)
    virtual_experiment = VirtualLinearModelExperiment("virtual_experiment_quadratic_model_meta.yaml")
    # write experiments raw data to file
    virtual_experiment.write_data_to_yaml("virtual_experiment_quadratic_model_data.yaml")


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
