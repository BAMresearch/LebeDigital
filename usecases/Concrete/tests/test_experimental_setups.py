import numpy as np

import usecases.Concrete.simulation_model as Simulation

import pytest

# import warnings
# from ffc.quadrature.deprecation \
#     import QuadratureRepresentationDeprecationWarning
# warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)

def simple_simulation(parameters,experiment):

    problem = Simulation.Models.LinearElasticity(experiment, parameters)

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
def old_test_experiemental_setup(dim, get_experiment):

    parameters = Simulation.Parameters() # using the current default values

    parameters['dim'] = dim
    parameters['mesh_density'] = 2
    parameters['log_level'] = 'WARNING'

    experiment = get_experiment(parameters)

    simple_simulation(parameters, experiment)


def test_new_setup():

    parameters = Simulation.Parameters() # using the current default values

    parameters['dim'] = 3
    parameters['mesh_density'] = 15
    parameters['log_level'] = 'WARNING'
    parameters['E'] = 3000
    parameters['nu'] = 0.2
    parameters['radius'] = 75
    parameters['height'] = 300
    parameters['bc_setting'] = 'fixed'

    experiment = Simulation.Setups.ConcreteCylinderExperiment(parameters)

    problem = Simulation.Models.LinearElasticity(experiment, parameters, pv_name='testing_cylinder')

    # data for time stepping
    dt = 1200  # 20 min step
    time = dt * 3  # total simulation time in s

    # set time step
    problem.set_timestep(dt)  # for time integration scheme

    # initialize time
    t = dt  # first time step time

    for displ in [-10,-50,-100]:
    #while t <= time:  # time


        problem.experiment.apply_displ_load(displ)
        # solve temp-hydration-mechanics
        problem.solve(t=t)  # solving this
        problem.pv_plot(t)

        # prepare next timestep
        t += dt