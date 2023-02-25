from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.computation_GWP_mix import computation_GWP_mix
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest


def test_computation_GWP_mix():

    constituents = {'cement' : {'content': 1 * ureg('kg/m^3'), 'GWP': 1 * ureg('kg_CO2_eq/kg')}}

    gwp_mix = computation_GWP_mix(constituents)

    assert gwp_mix == 1 * ureg('kg_CO2_eq/m^3')
