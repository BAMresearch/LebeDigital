import numpy as np


class AsymptoticModel:
    """
    Define a standard forward model for a specific type of test (in this case an asymptotic model that represents two data series MatB and MatC)

    Note:
        MatB(x) = data points for MatB sensor, following a + b*x in the virtual experiment (w/o noise)
        MatC(x) = data points for MatC sensor, following c + d*x + e*x² in the virtual experiment (w/o noise)
        x are the sensor positions in the interval [start_x, end_x].
        Different MatB and MatC have identical  sensor positions and sensor noise sepcifications.
        Model: Mat_B_C(x) = mA + mB / x, with Mat_B_C(x) = MatB(x)/MatC(x),
               mA and mB are model parameters to be calibrated.

    Attributes:
        x_Mat_sensors (np.array): Array of positions x where MatB and MatC values are measured
    """

    def __init__(self, x_Mat_sensors): #, mA, mB): #a, b, c, d, e
        """Construct an asymptotic forward model with sensor positions for MatB and MatC measurements"""
        self.x_Mat_sensors = x_Mat_sensors
        #self.a = a
        #self.b = b
        #self.c = c
        #self.d = d
        #self.e = e
#        self.mA = mA
#        self.mB = mB
#        self.x_MatC = x_MatC
#        self.x_Mat_B_C = x_Mat_B_C # needed?


    def __call__(self, parameters):
        """Evaluate the asymptotic model at predefined sensor positions for MatB and MatC values
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes (mA and mB)

        Returns:
            [MatB_x, MatC_x, Mat_B_C]: MatB, MatC (virtual experiment) and Mat_B_C (model) values at the locations defined in the constructor
        """
        self.check_parameters(parameters)
        # MatB_x = parameters['a'] * np.ones(len(self.x_Mat_sensors)) + parameters['b'] * self.x_Mat_sensors
        # MatC_x = parameters['c'] * np.ones(len(self.x_Mat_sensors)) + parameters['d'] * self.x_Mat_sensors + parameters['e'] * self.x_MatC * self.x_Mat_sensors
        #MatB_x = self.a * np.ones(len(self.x_Mat_sensors)) + self.b * self.x_Mat_sensors
        #MatC_x = self.c * np.ones(len(self.x_Mat_sensors)) + self.d * self.x_Mat_sensors + self.e * self.x_Mat_sensors * self.x_Mat_sensors

        Mat_B_C = parameters['mA'] * np.ones(len(self.x_Mat_sensors))  + parameters['mB'] / self.x_Mat_sensors 

#        print('DEBUG check_parameters in asymptotic_model.py:',parameters) # this is iterations

        return [Mat_B_C] #was also MatB_x, MatC_x,# rather not mA, mB

    @staticmethod
    def are_parameters_valid(parameters):
        """
        Check the parameters (all required are given and are also within the bounds)
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes
        Returns:
            true: if check is successful, other wise return false
        """
        try:
            AsymptoticModel.check_parameters()
        except Exception as exc:
            print(exc)
            return False
        return True

    @staticmethod
    def check_parameters(parameters):
        """
        Check the parameters (all required are given and are also within the bounds)
        Args:
            parameters(ModelParameters): a dictionary type of parameter list defined in bayes
        Raises:
            throws an exception when parameters are either not given or
            out of bounds
        """
        ## dont know, maybe also need a,b,c,d,e here. probably not.
        # could be.

        #assert 'a' in parameters.names
        #assert 'b' in parameters.names
        #assert 'c' in parameters.names
        #assert 'd' in parameters.names
        #assert 'e' in parameters.names
        assert 'mA' in parameters.names
        assert 'mB' in parameters.names
        #FINE # print('DEBUG param names',parameters.names)

        assert parameters['mA'] != 0 # mB doesnt work, assertion error.  mB not defined somewhere?!
        assert parameters['mB'] != 0 # should be in later
