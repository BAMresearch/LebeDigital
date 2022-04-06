import numpy as np

import usecases.Concrete.simulation_model as Simulation

import pytest

def simple_setup(p, displacement, sensor):
    parameters = Simulation.Parameters()  # using the current default values

    parameters['mesh_density'] = 4
    parameters['log_level'] = 'WARNING'
    parameters['bc_setting'] = 'free'
    parameters['dim'] = 2

    parameters = parameters + p

    experiment = Simulation.Setups.ConcreteCylinderExperiment(parameters)

    problem = Simulation.Models.LinearElasticity(experiment, parameters)
    problem.add_sensor(sensor)

    problem.experiment.apply_displ_load(displacement)

    problem.solve()  # solving this
    problem.pv_plot()

    # last measurement
    return problem.sensors[sensor.name].data[-1]



def test_poissions_ratio():
    p = Simulation.Parameters()  # using the current default values

    p['E'] = 2000
    p['nu'] = 0.499999999
    p['radius'] = 5
    p['height'] = 10
    displacement = -1

    sensor = Simulation.sensors.DisplacementSensor((p.radius*2,0))
    measured = simple_setup(p, displacement, sensor)

    print(measured)
    print(-p.nu*p.radius*2*displacement/p.height)
    #assert measured == pytest.approx(p.E*p.radius*2*displacement/p.height)


# testing the linear elastic response
def test_force_response():
    p = Simulation.Parameters()  # using the current default values

    p['E'] = 1023
    p['nu'] = 0.0
    p['radius'] = 6
    p['height'] = 12
    displacement = -3

    sensor = Simulation.sensors.ReactionForceSensor()
    measured = simple_setup(p, displacement, sensor)

    assert measured == pytest.approx(p.E*p.radius*2*displacement/p.height)


print('Moin')
test_poissions_ratio()
