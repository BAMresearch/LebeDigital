import pytest
from pint.testsuite.helpers import \
    assert_quantity_almost_equal as assert_approx

from lebedigital.demonstrator_scripts.computation_GWP_steel_per_volume import \
    computation_GWP_steel_per_volume
from lebedigital.unit_registry import ureg


def test_computation_GWP_per_part():
    # the values are chosen a priory and are just for testing purposes
    gwp_steel = 1 * ureg("kg_CO2_eq/kg")
    density_steel = 1 * ureg("kg/m^3")

    steel_gwp = computation_GWP_steel_per_volume(gwp_steel, density_steel)

    assert steel_gwp.magnitude == pytest.approx(1, 0.001)
    assert steel_gwp.units == ureg("kg_CO2_eq/m^3")
