import numpy as np

import fenics_concrete

import pytest


def test_cylinder_simulation():
    """Tesing the cylinder simulation
    
    This test is checking if the conda package is correctly installed"""

    parameters = fenics_concrete.Parameters()  # using the current default values

    parameters['E'] = 3000
    parameters['nu'] = 0.2
    parameters['radius'] = 75
    parameters['height'] = 300
    parameters['dim'] = 3

    parameters['log_level'] = 'WARNING'
    parameters['bc_setting'] = 'free'
    parameters['mesh_density'] = 6

    displacement = -3

    experiment = fenics_concrete.ConcreteCylinderExperiment(parameters)

    problem = fenics_concrete.LinearElasticity(experiment, parameters)
    sensor = fenics_concrete.sensors.ReactionForceSensorBottom()

    problem.add_sensor(sensor)

    problem.experiment.apply_displ_load(displacement)

    problem.solve()  # solving this

    # last measurement
    measured_force = problem.sensors[sensor.name].data[-1]

    # exact solution for free bc
    exact_solution = np.pi * parameters['radius'] ** 2 * parameters['E'] * displacement / parameters['height']

    #scaling factor due to an error resulting from the discretization
    rel_discretization_error = 1.012722

    assert measured_force == pytest.approx(rel_discretization_error*exact_solution)
