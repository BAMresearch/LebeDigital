from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.computation_GWP_per_part import computation_GWP_per_part
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest

def test_computation_GWP_per_part():

    # the values for GWP are chose a priory and are just for testing purposes
    gwp_concrete = 10 * ureg('kg_CO2_eq/m^3')
    width = 0.2 * ureg('m')
    height = 0.5 * ureg('m')
    length = 10 * ureg('m')
    gwp_steel = 1000 * ureg('kg_CO2_eq/m^3')
    n_steel = 2 * ureg('')
    diameter_steel = 10 * ureg('mm')

    beam_gwp = computation_GWP_per_part(gwp_concrete,gwp_steel,width,height,length,n_steel,diameter_steel)

    assert beam_gwp.magnitude == pytest.approx(11.570796326794897, 0.001)
    assert beam_gwp.units == ureg('kg_CO2_eq')

    concrete_gwp = computation_GWP_per_part(gwp_concrete,0 * ureg('kg_CO2_eq/m^3'),
                                            width,height,length,n_steel,diameter_steel)
    steel_gwp = computation_GWP_per_part(0 * ureg('kg_CO2_eq/m^3'),gwp_steel,
                                         width,height,length,n_steel,diameter_steel)

    assert beam_gwp.magnitude == pytest.approx((concrete_gwp + steel_gwp).magnitude)