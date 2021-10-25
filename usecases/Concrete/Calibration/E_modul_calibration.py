# -- PKM TUM -- atul.agrawal@tum.de -- #

from misc import load_experimental_data
import matplotlib.pyplot as plt

from probeye.definition.forward_model import ForwardModelBase
from probeye.definition.sensor import Sensor
from probeye.definition.inference_problem import InferenceProblem
from probeye.definition.noise_model import NormalNoiseModel
from probeye.inference.torch_.solver import run_pyro_solver
from probeye.postprocessing.sampling import create_trace_plot, create_pair_plot

import numpy as np
import matplotlib.pyplot as plt

# -- Loading the experiment
# for concrete E 60 - 85 Gpa/60 - 85 10^9N/m2 /60-85 KN/mm^2
skip_last = 145
skip_init = 330
exp_output = load_experimental_data('Wolf 8.2 Probe 1',skip_init,skip_last)
plt.plot(exp_output['stress'],exp_output['displacement']/exp_output['height'])  # checking loading data

# -- Set Numerical values
loc_E = 100
scale_E = 100 # this can be inferred too

# sigma prior parameters
low_sigma = 0
high_sigma = 0.005

# -- Define the forward solver
# TODO: to be replaced by PDE solver
class LinearModel(ForwardModelBase):
    def response(self,inp):
        x = inp['strain']
        E = inp['E']
        response = {}
        for os in self.output_sensors:
            e_hat = E  # rescaling
            response[os.name] = x*e_hat
        return response

# -- Set up inference problem
isensor = Sensor("strain")
osensor = Sensor("stress")
linear_model = LinearModel(['E'], [isensor], [osensor])

problem = InferenceProblem("Youngs Modulus Calibration")

problem.add_parameter('E', 'model',
                      tex="$E$",
                      info="Slope of the graph",
                      prior=('normal', {'loc': loc_E,
                                        'scale': scale_E}))

problem.add_parameter('sigma', 'noise',
                      tex=r"$\sigma$",
                      info="Std. dev, of 0-mean noise model",
                      prior=('uniform', {'low': low_sigma,
                                         'high': high_sigma}))

problem.add_noise_model(NormalNoiseModel(prms_def={'sigma': 'std'}, sensors=osensor))
problem.add_forward_model("LinearModel", linear_model)
problem.add_experiment(f'TestSeries_1', fwd_model_name="LinearModel",
                       sensor_values={isensor.name: exp_output['displacement']/exp_output['height'],
                                      osensor.name: exp_output['stress']})

problem.info()

# -- Solve inference problem
pos = run_pyro_solver(problem,n_steps=1000,n_initial_steps=100)

# -- Visualisation
create_trace_plot(pos, problem)

create_pair_plot(pos,problem)