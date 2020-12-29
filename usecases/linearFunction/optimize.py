import unittest
import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
from linear_model_error import *
from bayes.parameters import *
from bayes.multi_model_error import MultiModelError

class MyModelError(MultiModelError):
    "concatenate multiple 'standardized' model errors into a joint class and take care of the shared parameters"
    def __init__(self):
        # For the inference, we combine them and use a 'key' to distinguish
        # e.g. "A" from the one model to "A" from the other one.
        super().__init__()
        yaml_file_list = Path(Path(__file__).parents[0]).glob('LinearModelExperiment_*.yaml')
        for yaml_file in yaml_file_list:
            me = LinearModelError(str(yaml_file))
            parameter = me.get_parameter_dict()
            parameter.define("b")
            key = self.add(me, parameter)

        self.join(shared='b')
        self.set_latent('b')

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)

    my_model_error = MyModelError()
    start_vector = np.array([0.7])
    residual_vector = my_model_error([3])
    #print("residual", residual_vector)
    result = least_squares(my_model_error, start_vector)
    print(f"optimal set of parameters with b={result.x} and cost function {result.cost}.")