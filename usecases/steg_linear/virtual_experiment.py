import numpy as np
import yaml
from pathlib import Path


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
        self.sigma_noise_function = np.asarray(d['sigma_noise_function'])
        self.x_derivative = np.asarray(d['x_derivative'])
        self.sigma_noise_derivative = np.asarray(d['sigma_noise_derivative'])
        self.seed = np.asarray(d['seed'])

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
                    + self.b * self.x_function + self.c * np.square(self.x_function) \
                    + np.random.normal(0, self.sigma_noise_function, len(self.x_function))
                df_x = self.b * np.ones(len(self.x_derivative)) + 2. * self.c * self.x_derivative
                + np.random.normal(0, self.sigma_noise_derivative, len(self.x_derivative))

                data[index] = {
                    "a": float(a),
                    "f": f_x.tolist(),
                    "df": df_x.tolist()
                }
            yaml.dump(data, f, default_flow_style=None)


def main():
    metadata_files = Path(__file__).parent.glob("*_meta.yaml")
    for file in metadata_files:
        virtual_experiment = VirtualLinearModelExperiment(file)
        virtual_experiment.write_data_to_yaml(Path(str(file).replace('meta.yaml', 'data.yaml')))


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
