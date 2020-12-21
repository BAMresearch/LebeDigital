import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
from LinearModelError import *
from bayes.parameters import *

class MultiModelError:
    def __init__(self):
        self.prms = {}
        self.mes = {}
        self.n = 0
        self.joint_parameter_list = None
        pass

    def add(self, model_error, parameters, key=None):
        key = key or self.n
        assert key not in self.prms.keys()
        self.n += 1

        self.mes[key] = model_error
        self.prms[key] = parameters
        return key


    def __call__(self, parameter_vector):
        self.joint_parameter_list.update(parameter_vector)
        result = []
        for key, me in self.mes.items():
            prm = self.prms[key]
            result.append(me(prm))
        return np.concatenate(result)

    def join(self, shared=None):
        self.joint_parameter_list = JointParameterList(self.prms, shared)

    def set_latent(self, name, key=None):
        if self.joint_parameter_list is None:
            raise RuntimeError(
                "Error, you should first join all model errors to a joint list by calling join."
            )
        self.joint_parameter_list.set_latent(name, key)


class MyModelError:
    "concatenate multiple 'standardized' model errors into a joint class and take care of the shared parameters"
    def __init__(self):
        # For the inference, we combine them and use a 'key' to distinguish
        # e.g. "A" from the one model to "A" from the other one.
        self.me = MultiModelError()
        yaml_file_list = Path(Path(__file__).parents[0]).glob('LinearModelExperiment_*.yaml')
        for yaml_file in yaml_file_list:
            me = LinearModelError(str(yaml_file))
            parameter = me.get_parameter_dict()
            parameter.define("b")
            key = self.me.add(me, parameter)

        self.me.join(shared='b')
        self.me.set_latent('b')

    def __call__(self, parameter_vector):
        return self.me(parameter_vector)

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)

    my_model_error = MyModelError()
    start_vector = np.array([0.7])
    residual_vector = my_model_error([3])
    #print("residual", residual_vector)
    result = least_squares(my_model_error, start_vector)
    print(f"optimal set of parameters with b={result.x} and cost function {result.cost}.")