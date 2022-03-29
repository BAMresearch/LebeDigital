import dolfin as df
import numpy as np


class Sensors(dict):
    """
    Dict that also allows to access the parameter
        p["parameter"]
    via the matching attribute
        p.parameter
    to make access shorter
    """
    # TESTING a sensor dictionary
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        assert key in self
        self[key] = value

    def __setitem__(self, initial_key, value):
        # check if key exisits, if so, add a number
        i = 2
        key = initial_key
        if key in self:
            while key in self:
                key = initial_key + str(i)
                i += 1

        super().__setitem__(key, value)


# sensor template
class Sensor:
    def measure(self, u):
        raise NotImplementedError()

    @property
    def name(self):
        return self.__class__.__name__

    def data_max(self, value):
        if value > self.max:
            self.max = value


class DisplacementSensor(Sensor):
    def __init__(self, where):
        self.where = where
        self.data = []
        self.time = []

    def measure(self, problem, t=1.0):
        # get displacements
        self.data.append(problem.displacement(self.where))
        self.time.append(t)


class TemperatureSensor(Sensor):
    # temperature sensor in celsius
    def __init__(self, where):
        self.where = where
        self.data = []
        self.time = []

    def measure(self, problem, t=1.0):
        T = problem.temperature(self.where) - problem.p.zero_C
        self.data.append(T)
        self.time.append(t)


class MaxTemperatureSensor(Sensor):
    def __init__(self):
        self.data = [0.0]
        self.time = [0.0]
        self.max = 0.0

    def measure(self, problem, t=1.0):
        max_T = np.amax(problem.temperature.vector().get_local()) - problem.p.zero_C
        self.data.append(max_T)
        self.data_max(max_T)


class DOHSensor(Sensor):
    def __init__(self, where):
        self.where = where
        self.data = []
        self.time = []

    def measure(self, problem, t=1.0):
        # get DOH
        # TODO: problem with projected field onto linear mesh!?!
        alpha = problem.degree_of_hydration(self.where)
        self.data.append(alpha)
        self.time.append(t)


class MinDOHSensor(Sensor):
    def __init__(self):
        self.data = []
        self.time = []

    def measure(self, problem, t=1.0):
        # get min DOH
        min_DOH = np.amin(problem.q_degree_of_hydration.vector().get_local())
        self.data.append(min_DOH)
        self.time.append(t)


class MaxYieldSensor(Sensor):
    def __init__(self):
        self.data = [0.0]
        self.time = [0.0]
        self.max = 0.0

    def measure(self, problem, t=1.0):
        max_yield = np.amax(problem.q_yield.vector().get_local())
        self.data.append(max_yield)
        self.time.append(t)
        self.data_max(max_yield)
        
# TODO: add more sensor for other fields

