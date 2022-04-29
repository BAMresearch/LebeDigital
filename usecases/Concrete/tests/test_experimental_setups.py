import numpy as np

import usecases.Concrete.simulation_model as Simulation

import pytest


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("boundary_setting", ['free', 'fixed'])
def test_experiemental_setup(dim, boundary_setting):
    parameters = Simulation.Parameters()  # using the current default values

    parameters['mesh_density'] = 4
    parameters['log_level'] = 'WARNING'
    parameters['E'] = 3000
    parameters['nu'] = 0.2
    parameters['radius'] = 75
    parameters['height'] = 300
    parameters['bc_setting'] = boundary_setting
    parameters['dim'] = dim

    displacement = -parameters.height * 0.1

    experiment = Simulation.Setups.ConcreteCylinderExperiment(parameters)

    problem = Simulation.Models.LinearElasticity(experiment, parameters)

    problem.experiment.apply_displ_load(displacement)

    problem.solve()  # solving this
