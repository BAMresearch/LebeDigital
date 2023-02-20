import pytest
import fenics_concrete
from lebedigital.simulation.three_point_bending_beam import three_point_bending_example
from lebedigital.unit_registry import ureg

def test_three_point_bending_beam():
    # setup parameters:
    parameters = fenics_concrete.Parameters()  # using the current default values

    # input values for the material
    parameters['E'] = 30000 * ureg('N/mm^2')
    parameters['nu'] = 0.2 * ureg('')

    # run simulation
    stress_x = three_point_bending_example(parameters)

    assert stress_x.magnitude == pytest.approx(118.45516228258175)

    # checking for linearity
    multiplier = 2
    parameters['E'] = parameters['E'] * multiplier

    # run simulation
    stress_x_2 = three_point_bending_example(parameters)

    assert stress_x.magnitude*multiplier == pytest.approx(stress_x_2.magnitude)

