import numpy as np

import usecases.Concrete.simulation_model as Simulation

import pytest

# import warnings
# from ffc.quadrature.deprecation \
#     import QuadratureRepresentationDeprecationWarning
# warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)

def simple_simulation(parameters,experiment):

    problem = Simulation.Models.ConcreteThermoMechanical(experiment, parameters)

    # data for time stepping
    dt = 1200  # 20 min step
    time = dt * 3  # total simulation time in s

    # set time step
    problem.set_timestep(dt)  # for time integration scheme

    # initialize time
    t = dt  # first time step time

    while t <= time:  # time
        # solve temp-hydration-mechanics
        problem.solve(t=t)  # solving this

        # prepare next timestep
        t += dt


# testing the different experimental setup with options
# just checking that they run!

@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("get_experiment", [Simulation.Setups.ConcreteCubeExperiment,
                                            Simulation.Setups.MinimalCubeExperiment,
                                            Simulation.Setups.ConcreteColumnExperiment,
                                            Simulation.Setups.ConcreteBeamExperiment,
                                            ])
def test_experiemental_setup(dim, get_experiment):

    parameters = Simulation.Parameters() # using the current default values

    parameters['dim'] = dim
    parameters['mesh_density'] = 2
    parameters['log_level'] = 'WARNING'

    experiment = get_experiment(parameters)

    simple_simulation(parameters, experiment)

