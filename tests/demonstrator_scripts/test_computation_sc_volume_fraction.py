import pytest

from lebedigital.demonstrator_scripts.computation_sc_volume_fraction import computation_sc_volume_fraction
from lebedigital.unit_registry import ureg


def test_computation_sc_volume_fraction():
    # test the function with the following input
    # initialize input dictionary
    input_dic = {}

    # testing trivial case
    # densities
    input_dic["density_cem"] = 1000 * ureg("kg/m^3")
    input_dic["density_sub"] = 1000 * ureg("kg/m^3")

    # mass fractions
    input_dic["sc_mass_fraction"] = 0.452 * ureg("")
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    assert sc_volume_fraction == pytest.approx(input_dic["sc_mass_fraction"], 0.001)

    # mass fractions
    input_dic["sc_mass_fraction"] = 0.21754 * ureg("")
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    assert sc_volume_fraction == pytest.approx(input_dic["sc_mass_fraction"], 0.001)

    # mass fractions
    input_dic["sc_mass_fraction"] = 0.0 * ureg("")
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    assert sc_volume_fraction == pytest.approx(input_dic["sc_mass_fraction"], 0.001)

    # mass fractions
    input_dic["sc_mass_fraction"] = 1.0 * ureg("")
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    assert sc_volume_fraction == pytest.approx(input_dic["sc_mass_fraction"], 0.001)

    # now with different densities
    # densities
    input_dic["density_cem"] = 1000 * ureg("kg/m^3")
    input_dic["density_sub"] = 2000 * ureg("kg/m^3")

    # mass fractions
    input_dic["sc_mass_fraction"] = 0.5 * ureg("")
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    assert sc_volume_fraction == pytest.approx(1 / 3, 0.001)
