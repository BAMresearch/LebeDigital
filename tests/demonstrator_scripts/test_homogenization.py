import pytest
from lebedigital.simulation.concrete_homogenization import concrete_homogenization

def test_homogenization():
    # initialize dictionary
    parameters = {}

    # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
    # kg - m - s - N  - Pa - J

    # values are kind of made up but within the expected magnitude

    # paste data
    parameters['paste_E'] = 30e9  # Pa
    parameters['paste_nu'] = 0.2
    parameters['paste_C'] = 870  # J/kg Specific Heat Capacity
    parameters['paste_kappa'] = 1.8  # W/m/K Thermal conductivity
    parameters['paste_rho'] = 2400  # kg/m^3
    parameters['paste_fc'] = 30e6  # Pa
    parameters['paste_Q'] = 250000 # J/kg

    # aggregate data
    parameters['aggregates_E'] = 25e9  # Pa
    parameters['aggregates_nu'] = 0.3
    parameters['aggregates_C'] = 840  # J/kg Specific Heat Capacity
    parameters['aggregates_kappa'] = 0.8  # W/m/K Thermal conductivity
    parameters['aggregates_rho'] = 2600  # kg/m^3
    parameters['aggregates_vol_frac'] = 0.6

    results = concrete_homogenization(parameters)

    assert results['E'] == pytest.approx(27014932516.511917)
    assert results['nu'] == pytest.approx(0.26409495548961426)
    assert results['fc'] == pytest.approx(27652173.91304348)
    assert results['C'] == pytest.approx(2145600.0)
    assert results['rho'] == pytest.approx(2520.0)
    assert results['kappa'] == pytest.approx(1.152)
    assert results['Q'] == pytest.approx(240000000)