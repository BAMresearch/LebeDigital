import numpy as np
import yaml
from pathlib import Path

import bayes.correlation


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
        l_noise_function: correlation length of the noise for the function values
        x_derivative(np.array): positions of derivative sensors in the interval [0,1]
        l_noise_derivative: correlation length of the noise for the derivative values
        seed: seed of the random number generator
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
        self.l_noise_function = np.asarray(d['l_noise_function'])
        self.x_derivative = np.asarray(d['x_derivative'])
        self.sigma_noise_derivative = np.asarray(d['sigma_noise_derivative'])
        self.l_noise_derivative = np.asarray(d['l_noise_function'])
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
            np.random.seed(self.seed)

            for index, a in enumerate(self.all_a):
                if (self.l_noise_function > 0):
                    assert (len(self.x_function) < 201)
                    # interpolate the function on a fine grid with 1000 data points
                    # so that the noise distribution is not dependent on the discretization x
                    fine_grid = np.linspace(np.min(self.x_function), np.max(self.x_function), 1000)
                    noise_fine_grid = np.random.multivariate_normal(
                        np.zeros(len(fine_grid)),
                        bayes.correlation.squared_exponential(fine_grid,
                                                              self.l_noise_function) * self.sigma_noise_function ** 2, 1)[0]
                    noise_function = np.interp(self.x_function, fine_grid, noise_fine_grid)
                else:
                    noise_function = np.random.normal(0, self.sigma_noise_function, len(self.x_function))

                f_x = a * np.ones(len(self.x_function)) \
                    + self.b * self.x_function + self.c * np.square(self.x_function) + noise_function

                if (self.l_noise_derivative > 0):
                    assert (len(self.x_derivative) < 101)
                    # interpolate the function on a fine grid with 1000 data points
                    # so that the noise distribution is not dependent on the discretization x
                    fine_grid = np.linspace(np.min(self.x_derivative), np.max(self.x_derivative), 500)
                    noise_fine_grid = np.random.multivariate_normal(
                        np.zeros(len(fine_grid)),
                        bayes.correlation.squared_exponential(
                            fine_grid,
                            self.l_noise_derivative) * self.sigma_noise_derivative ** 2, 1)[0]
                    noise_derivative = np.interp(self.x_derivative, fine_grid, noise_fine_grid)
                else:
                    noise_derivative = np.random.normal(0, self.sigma_noise_derivative, len(self.x_derivative))

                df_x = self.b * np.ones(len(self.x_derivative)) + 2. * self.c * self.x_derivative
                + noise_derivative

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
