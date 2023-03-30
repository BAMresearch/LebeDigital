from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.computation_GWP_per_part import computation_GWP_per_part
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest

def test_computation_GWP_per_part():

    gwp = 10 * ureg('kg_CO2_eq/m^3')
    width = 0.2 * ureg('m')
    height = 0.5 * ureg('m')
    length = 10 * ureg('m')

    beam_gwp = computation_GWP_per_part(gwp,width,height,length)

    assert beam_gwp == gwp*width*height*length
