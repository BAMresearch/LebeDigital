import pytest
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.demonstrator_scripts.computation_specific_heat_capacity_paste import (
    computation_specific_heat_capacity_paste,
)
from lebedigital.unit_registry import ureg


def test_computation_GWP_per_part():
    # the values are chosen a priory and are just for testing purposes
    vol_frac_cement = 0.3 * ureg("")
    vol_frac_sub = 0.2 * ureg("")
    vol_frac_water = 0.5 * ureg("")
    shc_cement = 1 * ureg("J/kg/K")
    shc_sub = 2 * ureg("J/kg/K")
    shc_water = 3 * ureg("J/kg/K")
    density_cem = 1 * ureg("kg/m^3")
    density_sub = 1 * ureg("kg/m^3")
    density_water = 1 * ureg("kg/m^3")

    shc_paste = computation_specific_heat_capacity_paste(
        vol_frac_cement,
        vol_frac_sub,
        vol_frac_water,
        shc_cement,
        shc_sub,
        shc_water,
        density_cem,
        density_sub,
        density_water,
    )

    assert shc_paste.magnitude == pytest.approx(2.2, 0.001)
    assert shc_paste.units == ureg("J/kg/K")

    # test 2
    vol_frac_cement = 1 / 3 * ureg("")
    vol_frac_sub = 1 / 3 * ureg("")
    vol_frac_water = 1 / 3 * ureg("")
    shc_cement = 1 * ureg("J/kg/K")
    shc_sub = 1 * ureg("J/kg/K")
    shc_water = 1 * ureg("J/kg/K")
    density_cem = 1 * ureg("kg/m^3")
    density_sub = 1 * ureg("kg/m^3")
    density_water = 1 * ureg("kg/m^3")

    shc_paste = computation_specific_heat_capacity_paste(
        vol_frac_cement,
        vol_frac_sub,
        vol_frac_water,
        shc_cement,
        shc_sub,
        shc_water,
        density_cem,
        density_sub,
        density_water,
    )

    assert shc_paste.magnitude == pytest.approx(1, 0.001)
    assert shc_paste.units == ureg("J/kg/K")

    # test 3
    vol_frac_cement = 1 / 2 * ureg("")
    vol_frac_sub = 1 / 4 * ureg("")
    vol_frac_water = 1 / 4 * ureg("")
    shc_cement = 2 * ureg("J/kg/K")
    shc_sub = 1 * ureg("J/kg/K")
    shc_water = 1 * ureg("J/kg/K")
    density_cem = 1 * ureg("kg/m^3")
    density_sub = 1 * ureg("kg/m^3")
    density_water = 1 * ureg("kg/m^3")

    shc_paste = computation_specific_heat_capacity_paste(
        vol_frac_cement,
        vol_frac_sub,
        vol_frac_water,
        shc_cement,
        shc_sub,
        shc_water,
        density_cem,
        density_sub,
        density_water,
    )

    assert shc_paste.magnitude == pytest.approx(1.5, 0.001)
    assert shc_paste.units == ureg("J/kg/K")
