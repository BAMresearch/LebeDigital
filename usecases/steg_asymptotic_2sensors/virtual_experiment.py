import numpy as np
import yaml
from pathlib import Path


class VirtualAsymptoticModelExperiment:
    """create virtual test data for a MatB and MatC

    Note:
        The model is comprised of 2 data series:
       MatB = a + b*x and MatC = c + d*x + e*x². a,b,c,d,e are used to generate the model data.
      The target property for the model to be calibrated is MatB/MatC = mA * mB / x, where mA and mB are model parameters.

    Attributes:
        a, b, c, d, e (floats) : coefficient of the virtual data series.
        mA (float): coefficient of the model
        mB (float): coefficient of the model
        x_Mat_sensors(np.array): positions of MatB and MatC sensors in the interval [x_start, x_end]

    """

    def __init__(self, virtual_experiment_metadata_yaml):
        """Create virtual experiment
        Args:
            virtual_experiment_metadata_yaml: meta data file to generate the data from
        """
        with open(virtual_experiment_metadata_yaml, "r") as f:
            dat = yaml.load(f, Loader=yaml.FullLoader)
        self.a = dat['a']
        self.b = dat['b']
        self.c = dat['c']
        self.d = dat['d']
#        self.e = dat['e']
      #  self.mA = dat['mA'] # maybe not
       # self.mB = dat['mB'] # maybe not
        self.x_Mat_sensors = np.asarray(dat['x_Mat_sensors'])
        self.sigma_noise_Mat = np.asarray(dat['sigma_noise_Mat']) # only one sensor with one noise parameter for both measurements, MatB and MatC.
        self.seed = np.asarray(dat['seed'])

    def write_data_to_yaml(self, virtual_experimental_data_file):
        """Write virtual sensor data to yaml file

        Args:
            virtual_experimental_data_file(str): file location

        Returns:
            writes the experimental data to the yaml file
        """
        with open(virtual_experimental_data_file, "w") as f:
            data = {}
            MatB_x = self.a * np.ones(len(self.x_Mat_sensors)) \
                    + self.b * self.x_Mat_sensors \
                    + np.random.normal(0, self.sigma_noise_Mat, len(self.x_Mat_sensors))
            MatC_x = self.c * np.ones(len(self.x_Mat_sensors)) \
                    + self.d*  self.x_Mat_sensors  + self.c * np.square(self.x_Mat_sensors) \
                    + np.random.normal(0, self.sigma_noise_Mat, len(self.x_Mat_sensors))
            Mat_B_C_x = MatB_x / MatC_x 
            data= {
                    "MatB": MatB_x.tolist(),
                    "MatC": MatC_x.tolist(),
                    "Mat_B_C": Mat_B_C_x.tolist()  
                }
            yaml.dump(data, f, default_flow_style=None)


def main():
    metadata_files = Path(__file__).parent.glob("*_meta.yaml")
    for file in metadata_files:
        virtual_experiment = VirtualAsymptoticModelExperiment(file)
        virtual_experiment.write_data_to_yaml(Path(str(file).replace('meta.yaml', 'data.yaml')))


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    main()
