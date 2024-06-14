from fenics import *
import numpy as np
import matplotlib.pyplot as plt
#import concrete_experiment as concrete_experiment
import fenics_concrete

# import probeye
from probeye.definition.inference_problem import InferenceProblem
from probeye.definition.forward_model import ForwardModelBase
from probeye.definition.sensor import Sensor
from probeye.definition.likelihood_model import GaussianLikelihoodModel
from probeye.inference.scipy.solver import ScipySolver


# ============================================================================ #
#                           Define the Forward Model                           #
# ==========================================
class HydrationHeatModelStep(ForwardModelBase):
    def definition(self):
        self.parameters = ['eta', 'B1', 'B2', "E_act"]
        # irgendeine liste....
        self.input_sensors = [Sensor("T"),
                              Sensor("dt"),
                              Sensor("time"),
                              # Sensor("E_act"),
                              Sensor("Q_pot"),
                              Sensor("T_ref"),
                              Sensor("alpha_max"),
                              Sensor("time")]
        self.output_sensors = [Sensor('heat')]

    def response(self, inp: dict) -> dict:
        # this method *must* be provided by the user
        T = inp["T"]
        dt = inp["dt"]
        time_list = inp["time"]
        parameter = {}
        parameter['B1'] = inp["B1"]
        parameter['B2'] = inp["B2"]
        parameter['eta'] = inp["eta"]
        parameter['alpha_max'] = inp["alpha_max"]
        parameter['E_act'] = inp["E_act"]
        parameter['T_ref'] = inp["T_ref"]
        parameter['Q_pot'] = inp["Q_pot"]

        # initiate material problem
        material_problem = fenics_concrete.ConcreteThermoMechanical()
        # get the respective function
        hydration_fkt = material_problem.get_heat_of_hydration_ftk()

        heat_list, dummy = hydration_fkt(T, time_list, dt, parameter)
        return {'heat': heat_list}


# ------------------------------------------
# START PROBLEM DESCRIPTION!!!!!!!
# -------------------------------------------
# read data
time_data = []
heat_data = []

T_datasets = []

# extract data from csv file
with open('cost_action_hydration_data.csv') as f:
    for i, line in enumerate(f):
        if i == 0:
            split_line = line.split(',')
            for j in range(0, len(split_line), 2):
                degree = split_line[j].split('_')[0]
                T_datasets.append(float(degree.strip()))
                time_data.append([])
                heat_data.append([])
        if i > 1:
            split_line = line.split(',')
            for j in range(len(T_datasets)):
                print(i, j, split_line[j * 2], split_line[j * 2 + 1])
                if split_line[j * 2].strip() != '':
                    time_data[j].append(float(split_line[j * 2].strip()) * 60 * 60)  # convert to seconds
                    heat_data[j].append(float(split_line[j * 2 + 1].strip()))

# sort data!!!
for i in range(len(heat_data)):
    zipped_lists = zip(time_data[i], heat_data[i])
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    time_data[i], heat_data[i] = [list(tuple) for tuple in tuples]

# ============================================================================ #
#                              Set numeric values                              #
# ============================================================================ #

problem = InferenceProblem("Linear regression with normal additive error")

problem.add_parameter(
    "eta",
    "model",
    tex=r"$\eta$",
    info="Some parameter, but important",
    prior=("normal", {"loc": 5.5, "scale": 1}),
)

problem.add_parameter(
    "B1",
    "model",
    tex=r"$B_1$",
    info="Some other parameter, but important",
    prior=("normal", {"loc": 0.00029, "scale": 0.001}),
    # prior=("uniform", {"low": 0.0, "high": 0.1}),
)

problem.add_parameter(
    "B2",
    "model",
    tex=r"$B_2$",
    info="Some other parameter, but important",
    prior=("normal", {"loc": 0.0024, "scale": 0.001}),
    # prior=("uniform", {"low": 0.0, "high": 1.0}),
    # prior=("uniform", {"low": 0.0, "high": 1.0}),
)

problem.add_parameter(
    "E_act",
    "model",
    tex=r"$E_act$",
    info="Some other parameter, but important",
    prior=("normal", {"loc": 47002, "scale": 10000}),
    # prior=("uniform", {"low": 0.0, "high": 1.0}),
    # prior=("uniform", {"low": 0.0, "high": 1.0}),
)

problem.add_parameter(
    "sigma",
    "likelihood",
    tex=r"$\sigma",
    info="Some parameter, but important",
    # prior=("uniform", {"low": 0.001, "high": 1}),
    const=0.01
)

hydration_heat_model = HydrationHeatModelStep()
problem.add_forward_model("HydrationHeatModel", hydration_heat_model)

# add the experimental data

for i, T in enumerate(T_datasets):
    problem.add_experiment(
        f"TestSeries_{i}",
        fwd_model_name="HydrationHeatModel",
        sensor_values={
            'time': time_data[i],
            'heat': heat_data[i],
            'alpha_max': 0.85,
            # 'E_act': 47002,   # activation energy in Jmol^-1
            # 'E_act': 42,   # dummy value for T = T_ref
            'T_ref': 25,  # reference temperature in degree celsius
            'Q_pot': 450e3,  # potential heat per weight of binder in J/kg
            'T': T,
            'dt': 300,
        },
    )

# add the noise model to the problem
problem.add_likelihood_model(
    GaussianLikelihoodModel(
        prms_def={"sigma": "std_model"}, sensors=[hydration_heat_model.output_sensors[0]]
    )
)

# give problem overview
problem.info()

# solve the thing!!!
scipy_solver = ScipySolver(problem)
inference_data = scipy_solver.run_max_likelihood(solver_options={"maxiter": 1000})

time_list = []
heat_list = []

# generate a time list for plotting
# get max time
tmax = 0
for i in range(len(time_data)):
    if time_data[i][-1] > tmax:
        tmax = time_data[i][-1]

dt = problem.experiments[f'TestSeries_{i}']['sensor_values']['dt']
plot_time_list = np.arange(0, tmax, dt)

for i, T in enumerate(T_datasets):
    time_list.append([])
    heat_list.append([])

    # check results
    vars = problem.experiments[f'TestSeries_{i}']['sensor_values']
    # set required parameter
    parameter = {}  # using the current default values
    parameter['B1'] = inference_data.x[1]  # in 1/s (le 0, smaller 0.1)
    parameter['B2'] = inference_data.x[2]  # - (le 0, smaller 1)
    parameter['eta'] = inference_data.x[0]  # something about diffusion  (should be larger 0)
    parameter['alpha_max'] = vars[
        'alpha_max']  # also possible to approximate based on equation with w/c (larger 0 and max 1)
    parameter['E_act'] = inference_data.x[3]  # vars['E_act']   # activation energy in Jmol^-1 (no relevant limits)
    parameter['T_ref'] = vars['T_ref']  # reference temperature in degree celsius
    parameter['Q_pot'] = vars['Q_pot']  # potential heat per weight of binder in J/kg
    dt = vars['dt']
    time_list[i] = plot_time_list

    # initiate material problem
    material_problem = fenics_concrete.ConcreteThermoMechanical()
    # get the respective function
    hydration_fkt = material_problem.get_heat_of_hydration_ftk()
    heat_list[i], dummy = hydration_fkt(T, time_list[i], dt, parameter)

    plt.plot(time_list[i], heat_list[i], color='black')
    plt.plot(time_data[i], heat_data[i], color='red', linestyle='dashed')

plt.show()