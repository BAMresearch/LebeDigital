import numpy as np
from asymptotic_model import AsymptoticModel
from bayes.parameters import ModelParameters
from collections import OrderedDict


class AsymptoticModelError:
    """
    Define a standard model error for a specific type of test (in this case a asymptotic model)

    Note:
        f(x) = mA + mB/x , where mA and mB are model parameters
        to be calibrated and x are the sensor positions in the interval [x_start, x_end].
        MatB and MatC have identical sensor positions and noise specs.
        MatB = a + bx  
        MatC = c + dx + ex²
        and Mat_B_C = MatB / MatC
        f(x): Mat_B_C = mA + mB/x

    Attributes:
        x_Mat_sensors (np.array): Array of positions x where Mat values are measured
        data_MatB (np.array): Array of sensor readings at x_Mat_sensors
        data_MatC (np.array): Array of sensor readings at x_Mat_sensors
        data_Mat_B_C  : (np.array): Array of evaluated MatB / MatC data
        asymptotic_model (AsymptoticModel): asymptotic forward model (mA+mB/x) , where 'mA' and 'mB' are model parameters
    """

    def __init__(self, x_Mat_sensors, data_MatB, data_MatC, data_Mat_B_C):
        """Create a asymptotic model error

        Args:
            x_Mat_sensors (np.array): Array of positions x where function values are measured
            data_MatB (np.array): Array of sensor readings at x_MatB
            data_MatC (np.array): Array of sensor reading at x_MatC
            data_Mat_B_C (np.array): Array of MatB / MatC evaluation
        """
        self.x_Mat_sensors = x_Mat_sensors
        self.data_MatB = data_MatB
        self.data_MatC = data_MatC
        self.data_Mat_B_C = data_Mat_B_C
#        self.a = a # rather not
#        self.b = b # rather not
#        self.c = c # rather not
#        self.d = d # rather not
#        self.e = e # rather not
        self.asymptotic_model = AsymptoticModel(self.x_Mat_sensors)

    def __call__(self, parameters):
        """Evaluate the model error for a specific experiment

        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (a and b)

        Returns:
            concatenated model error (np.array) of all individual model errors
        """
        # [MatB, MatC, Mat_B_C] = self.asymptotic_model(parameters)
        [Mat_B_C] = self.asymptotic_model(parameters)
#        return np.concatenate((f-self.data_MatB, df-self.data_MatC))
        return np.concatenate((Mat_B_C - self.data_MatB / self.data_MatC ))

    def evaluate(self, parameters):
        """Evaluate the model error for a specific experiment

        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (a and b)

        Returns:
            return an OrderDict of the individual model errors (e.g. per sensor)
        """
        [Mat_B_C] = self.asymptotic_model(parameters)
        return OrderedDict([('Mat_B_C', Mat_B_C-self.data_MatB/ self.data_MatC)])

    def get_parameter_dict(self):
        """Create a parameter list initialized with parameters given in the experimental data file

        Returns: parameter list
        """
        prm = ModelParameters()
        prm.define("a", self.a) # old
        # prm.define("mA",self.mA) # maybe?
        return prm
