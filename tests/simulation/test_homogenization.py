import pytest
from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

def test_homogenization():
    # initialize dictionary
    parameters = {}

    # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
    # kg - m - s - N  - Pa - J

    # values are kind of made up but within the expected magnitude

    # paste data
    parameters['paste_E'] = 30e9 * ureg('Pa')
    parameters['paste_nu'] = 0.2
    parameters['paste_C'] = 870 * ureg('J/kg')  # Specific Heat Capacity
    parameters['paste_kappa'] = 1.8 * ureg('W/m/K')  # Thermal conductivity
    parameters['paste_rho'] = 2400 * ureg('kg/m^3')
    parameters['paste_fc'] = 30e6 * ureg('Pa')
    parameters['paste_Q'] = 250000 * ureg('J/kg')

    # aggregate data
    parameters['aggregates_E'] = 25e9 * ureg('Pa')
    parameters['aggregates_nu'] = 0.3
    parameters['aggregates_C'] = 840 * ureg('J/kg')  # Specific Heat Capacity
    parameters['aggregates_kappa'] = 0.8 * ureg('W/m/K')  # Thermal conductivity
    parameters['aggregates_rho'] = 2600 * ureg('kg/m^3')
    parameters['aggregates_vol_frac'] = 0.6

    results = concrete_homogenization(parameters)

    assert_approx(results['E'], 27014932516.511917 * ureg('Pa'), rtol=0.001)
    assert results['nu'] == pytest.approx(0.26409495548961426)
    assert_approx(results['fc'], 27652173.91304348 * ureg('Pa'), rtol=0.001)
    assert_approx(results['C'], 2145600.0 * ureg('J/m^3'), rtol=0.001)
    assert_approx(results['rho'], 2520.0* ureg('kg/m^3'), rtol=0.001)
    assert_approx(results['kappa'], 1.152 * ureg('W/m/K'), rtol=0.001)
    assert_approx(results['Q'], 240000000 * ureg('J/m^3'), rtol=0.001)