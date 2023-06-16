from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.demonstrator_scripts.computation_volume_content import computation_ratios, computation_volume_content
from lebedigital.unit_registry import ureg


def test_computation_volume_content():
    """
    This is the test for the computation_volume_content function.
    It tests all outputs and performs a number of sanity checks
    """

    # TEST No 1, check intended use in general
    input1 = {}
    # densities
    input1["density_cem"] = 1.44 * ureg("g/cm^3")
    input1["density_sub"] = 840 * ureg("kg/m^3")
    input1["density_water"] = 977 * ureg("kg/m^3")
    input1["density_plasticizer"] = 0.98 * ureg("g/cm^3")
    input1["density_aggregates"] = 1.5 * ureg("t/m^3")

    input1["wb_mass_ratio"] = 0.4 * ureg("dimensionless")
    input1["sc_mass_fraction"] = 0.4 * ureg("dimensionless")
    input1["aggregates_volume_fraction"] = 0.65 * ureg("dimensionless")
    input1["plasticizer_volume_content"] = 10 * ureg("kg/m^3")

    output1 = computation_volume_content(input1)

    # check all values
    assert_approx(output1["water_vol_fraction"], 0.09983100608664518 * ureg("dimensionless"), rtol=1e-6)
    assert_approx(output1["sub_vol_fraction"], 0.127981287 * ureg("dimensionless"), rtol=1e-6)
    assert_approx(output1["cem_vol_fraction"], 0.111983626 * ureg("dimensionless"), rtol=1e-6)
    assert_approx(output1["cem_vol_content"], 161.256421 * ureg("kg/m^3"), rtol=1e-6)
    assert_approx(output1["sub_vol_content"], 107.504281 * ureg("kg/m^3"), rtol=1e-6)
    assert_approx(output1["water_vol_content"], 97.5348929 * ureg("kg/m^3"), rtol=1e-6)
    assert_approx(output1["aggregates_vol_content"], 975.0 * ureg("kg/m^3"), rtol=1e-6)
    assert_approx(output1["density_paste"], 1075.13027 * ureg("kg/m^3"), rtol=1e-6)
    assert_approx(output1["sc_volume_fraction"], 0.533333333 * ureg(""), rtol=1e-6)

    # TEST No 2, sanity check density paste
    input2 = {}
    # densities
    input2["density_cem"] = 1000 * ureg("kg/m^3")
    input2["density_sub"] = 1000 * ureg("kg/m^3")
    input2["density_water"] = 1000 * ureg("kg/m^3")
    input2["density_plasticizer"] = 1000 * ureg("kg/m^3")
    input2["density_aggregates"] = 1000 * ureg("kg/m^3")

    input2["wb_mass_ratio"] = 0.4 * ureg("dimensionless")
    input2["sc_mass_fraction"] = 0.4 * ureg("dimensionless")
    input2["aggregates_volume_fraction"] = 0.65 * ureg("dimensionless")
    input2["plasticizer_content"] = 10 * ureg("kg/m^3")

    output2 = computation_volume_content(input2)

    # check paste density
    assert_approx(output2["density_paste"], 1000 * ureg("kg/m^3"), rtol=1e-6)

    # TEST No 3, sanity check w/b ratio of 1, compare weight with equal density
    input3 = {}
    # densities
    input3["density_cem"] = 1000 * ureg("kg/m^3")
    input3["density_sub"] = 840 * ureg("kg/m^3")
    input3["density_water"] = 1000 * ureg("kg/m^3")
    input3["density_plasticizer"] = 0.98 * ureg("g/cm^3")
    input3["density_aggregates"] = 1.5 * ureg("t/m^3")

    input3["wb_mass_ratio"] = 1
    input3["sc_mass_fraction"] = 0
    input3["aggregates_volume_fraction"] = 0

    output3 = computation_volume_content(input3)

    # check all values
    assert_approx(output3["water_vol_fraction"], output3["cem_vol_fraction"], rtol=1e-6)


def test_computation_slag_ratio():
    """
    testing minimal option of computation of the ratios
    """

    # compute ratios and fractions based on content, minimal input

    # TEST No 1, check the output of changes to sc_volume_fraction
    # reference
    input1 = {}
    # densities
    input1["density_cem"] = 1.44 * ureg("g/cm^3")
    input1["density_sub"] = 840 * ureg("kg/m^3")
    input1["density_water"] = 977 * ureg("kg/m^3")
    input1["density_plasticizer"] = 0.98 * ureg("g/cm^3")
    input1["density_aggregates"] = 1.5 * ureg("t/m^3")

    input1["wb_mass_ratio"] = 0.4 * ureg("dimensionless")
    input1["sc_mass_fraction"] = 0.4 * ureg("dimensionless")
    input1["aggregates_volume_fraction"] = 0.65 * ureg("dimensionless")
    input1["plasticizer_volume_content"] = 10 * ureg("kg/m^3")

    output1 = computation_volume_content(input1)

    # checking extremes
    sc_vol_fraction = output1["sub_vol_fraction"] + output1["cem_vol_fraction"]

    input2 = input1
    input2["sc_mass_fraction"] = 0.0 * ureg("dimensionless")
    output2 = computation_volume_content(input2)

    assert_approx(output2["sub_vol_fraction"], 0.0, rtol=0.05)
    assert_approx(output2["cem_vol_fraction"], sc_vol_fraction, rtol=0.1)

    input3 = input1
    input3["sc_mass_fraction"] = 1.0 * ureg("dimensionless")

    output3 = computation_volume_content(input3)

    assert_approx(output3["sub_vol_fraction"], sc_vol_fraction, rtol=0.1)
    assert_approx(output3["cem_vol_fraction"], 0.0, rtol=0.05)


def test_computation_ratio():
    """
    testing minimal option of computation of the ratios
    """

    # compute ratios and fractions based on content, minimal input
    input1 = {}
    # densities
    input1["density_cem"] = 1543 * ureg("kg/m^3")
    input1["density_water"] = 988 * ureg("kg/m^3")

    input1["cem_vol_content"] = 1 * ureg("kg/m^3")
    input1["water_vol_content"] = 1 * ureg("kg/m^3")

    output1 = computation_ratios(input1)

    assert_approx(output1["wb_mass_ratio"], 1, rtol=1e-6)
    assert_approx(output1["sc_volume_fraction"], 0, rtol=1e-6)
    assert_approx(output1["aggregates_volume_fraction"], 0, rtol=1e-6)


def test_computation_ratio_based_on_computation_volume_content():
    """
    testing the computation of the rations, based on the output of the volume computation
    """

    # compute volume content
    input1 = {}
    # densities
    input1["density_cem"] = 1543 * ureg("kg/m^3")
    input1["density_sub"] = 900 * ureg("kg/m^3")
    input1["density_water"] = 988 * ureg("kg/m^3")
    input1["density_plasticizer"] = 957 * ureg("kg/m^3")
    input1["density_aggregates"] = 1.46 * ureg("t/m^3")

    input1["wb_mass_ratio"] = 0.765
    input1["sc_mass_fraction"] = 0.452
    input1["aggregates_volume_fraction"] = 0.578
    input1["plasticizer_content"] = 12 * ureg("kg/m^3")

    output1 = computation_volume_content(input1)

    # compute ratios and fractions based on content
    input2 = {}
    # densities
    input2["density_cem"] = input1["density_cem"]
    input2["density_sub"] = input1["density_sub"]
    input2["density_water"] = input1["density_water"]
    input2["density_plasticizer"] = input1["density_plasticizer"]
    input2["density_aggregates"] = input1["density_aggregates"]

    input2["plasticizer_content"] = input1["plasticizer_content"]
    input2["cem_vol_content"] = output1["cem_vol_content"]
    input2["sub_vol_content"] = output1["sub_vol_content"]
    input2["water_vol_content"] = output1["water_vol_content"]
    input2["aggregates_vol_content"] = output1["aggregates_vol_content"]

    output2 = computation_ratios(input2)

    assert_approx(input1["wb_mass_ratio"], output2["wb_mass_ratio"], rtol=1e-6)
    assert_approx(output1["sc_volume_fraction"], output2["sc_volume_fraction"], rtol=1e-6)
    assert_approx(input1["aggregates_volume_fraction"], output2["aggregates_volume_fraction"], rtol=1e-6)
