# -- PKM TUM -- atul.agrawal@tum.de -- #
# local imports
from usecases.Concrete.Prediction.three_point_bending_example import three_point_bending_example
from misc import load_experimental_data, PosteriorPredictive
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from probeye.definition.forward_model import ForwardModelBase
from probeye.definition.sensor import Sensor
from probeye.definition.inference_problem import InferenceProblem
from probeye.definition.noise_model import NormalNoiseModel
from probeye.inference.torch_.solver import PyroSolver
from probeye.postprocessing.sampling import create_trace_plot, create_pair_plot


# -- Loading the experiment
# for concrete E 60 - 85 Gpa/60 - 85 10^9N/m2 /60-85 KN/mm^2
skip_last = 145
skip_init = 330
experiment_name = 'Wolf 8.2 Probe 1'
path = '../knowledgeGraph/emodul/E-modul-processed-data/processeddata/processed_Wolf_8_2_Probe_1.csv'
exp_output = load_experimental_data(experiment_name,skip_init,skip_last, path=path)
plt.plot(exp_output['stress'],exp_output['displacement']/exp_output['height'])  # checking loading data

# -- Set Numerical values
# "uninformed" Normal prior for the E modulus
loc_E = 100 # kN/mm2
scale_E = 100 # this can be inferred too

# Uniform prior for experimental noise/model discrepancy
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
solver = PyroSolver(problem)
pos = solver.run_mcmc(n_steps=500,n_initial_steps=100)
mcmc_pyro = solver.raw_results
mcmc_pyro.summary()

# -- Storing in KG with RDF
# TODO

# -- Visualisation
create_trace_plot(pos, problem)

create_pair_plot(pos,problem,focus_on_posterior=True)

# -- Posterior Predictive
# ---- query KG
# TODO

# ---- posterior predictive values
nu = 0.2
E_pos = 1000*solver.raw_results.get_samples()['E'].numpy() # N/mm2
three_point = three_point_bending_example()
pos_pred = PosteriorPredictive(three_point.run, known_input_solver=nu,parameter=E_pos)
mean, sd = pos_pred.get_stats(samples=20)
# ---- visualize posterior predictive
posterior_pred_samples = pos_pred._samples
np.save('./post_pred.npy',posterior_pred_samples)
sns.kdeplot(posterior_pred_samples)
plt.xlabel('stress in x-direction')


