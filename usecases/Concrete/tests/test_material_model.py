import numpy as np

import fenics_concrete

import pytest

def simple_setup(p, displacement, sensor):
    parameters = fenics_concrete.Parameters()  # using the current default values

    parameters['log_level'] = 'WARNING'
    parameters['bc_setting'] = 'free'
    parameters['mesh_density'] = 6

    parameters = parameters + p

    experiment = fenics_concrete.ConcreteCylinderExperiment(parameters)

    problem = fenics_concrete.LinearElasticity(experiment, parameters)
    problem.add_sensor(sensor)

    problem.experiment.apply_displ_load(displacement)

    problem.solve()  # solving this

    # last measurement
    return problem.sensors[sensor.name].data[-1]

# testing the linear elastic response
def test_force_response_2D():
    p = fenics_concrete.Parameters()  # using the current default values

    p['E'] = 1023
    p['nu'] = 0.0
    p['radius'] = 6
    p['height'] = 12
    displacement = -3
    p['dim'] = 2

    sensor = fenics_concrete.sensors.ReactionForceSensor()
    measured = simple_setup(p, displacement, sensor)

    assert measured == pytest.approx(p.E*p.radius*2*displacement/p.height)

def test_force_response_3D():
    p = fenics_concrete.Parameters()  # using the current default values

    p['E'] = 1023
    p['nu'] = 0.0
    p['radius'] = 6
    p['height'] = 12
    displacement = -3
    p['dim'] = 3

    sensor = fenics_concrete.sensors.ReactionForceSensor()
    measured = simple_setup(p, displacement, sensor)

    # due to meshing errors, only aprroximate results to be expected. within 1% is good enough
    assert measured == pytest.approx(p.E*np.pi*p.radius**2*displacement/p.height, 0.01)
