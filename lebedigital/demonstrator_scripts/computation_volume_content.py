import pytest

from lebedigital.unit_registry import ureg


def computation_volume_content(input_dic):
    """
    This is the function to compute volume contents based on mass or volume ratios

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Volume fractions refers to the fraction of the mix volume a phase is occupying.
    Volume content is the mass of the phase within a 1m^3.

    Parameters
    ----------
    input_dic : dic / pint units required converted to matching units
        required:
            - 'density_cem'
            - 'density_water'
            - 'wb_mass_ratio' : water to binder ratio water/(cement + substitute (e.g. slag or cemII))
            - 'density_aggregates'
            - 'aggregates_volume_fraction'
            - 'density_sub' : density of the substitute (slag or cemII)
            - 'sc_mass_fraction' : mass fraction of substitute to base cement
        optional (otherwise assumed to be zero)
             - 'plasticizer_volume_content'
             - 'density_plasticizer'
             - 'sc_volume_fraction' : volume fraction of substitute to base cement

    Returns
    -------
    output : dic / pint unit given
             - 'water_vol_fraction'
             - 'sub_vol_fraction'
             - 'cem_vol_fraction'
             - 'cem_vol_content'
             - 'sub_vol_content'
             - 'water_vol_content'
             - 'aggregates_vol_content'
             - 'density_paste'
    """

    # initialize output dictionary
    output = {}

    # optional parameter
    if "plasticizer_volume_content" not in input_dic.keys():
        input_dic["plasticizer_volume_content"] = 0 * ureg("kg/m^3")
        input_dic["density_plasticizer"] = 42 * ureg("kg/m^3")  # dummy value

    # converting to correct pint units / automatic check for pint input / check for required input values
    input_dic["density_sub"].ito("kg/m^3")
    input_dic["density_cem"].ito("kg/m^3")
    input_dic["density_water"].ito("kg/m^3")
    input_dic["density_plasticizer"].ito("kg/m^3")
    input_dic["density_aggregates"].ito("kg/m^3")

    # compute volume fraction of substitute to cement
    if input_dic["sc_mass_fraction"] == 0:
        output["sc_volume_fraction"] = 0
    else:
        mass_sub = 1 * ureg("kg") / (1 + (1 - input_dic["sc_mass_fraction"]) / input_dic["sc_mass_fraction"])
        mass_cem = 1 * ureg("kg") - mass_sub
        vol_sub = mass_sub / input_dic["density_sub"]
        vol_cem = mass_cem / input_dic["density_cem"]
        output["sc_volume_fraction"] = vol_sub / (vol_sub + vol_cem)

    # compute binder density

    density_binder = (
        output["sc_volume_fraction"] * input_dic["density_sub"]
        + (1 - output["sc_volume_fraction"]) * input_dic["density_cem"]
    )

    # compute volume ratio of water to binder (cement and slag)
    water_vol_fraction_to_binder = (
        input_dic["wb_mass_ratio"]
        * density_binder
        / (input_dic["density_water"] + input_dic["wb_mass_ratio"] * density_binder)
    )

    # volume fractions
    ## plasticizer
    output["plasticizer_vol_fraction"] = input_dic["plasticizer_volume_content"] / input_dic["density_plasticizer"]
    ## water
    output["water_vol_fraction"] = (
        1 - input_dic["aggregates_volume_fraction"]
    ) * water_vol_fraction_to_binder - output["plasticizer_vol_fraction"]
    ## binder
    binder_vol_fraction = (1 - input_dic["aggregates_volume_fraction"]) * (1 - water_vol_fraction_to_binder)

    output["sub_vol_fraction"] = binder_vol_fraction * output["sc_volume_fraction"]
    output["cem_vol_fraction"] = binder_vol_fraction * (1 - output["sc_volume_fraction"])

    # sanity check, the volume fractions add up to 1
    assert 1 == pytest.approx(
        (
            output["plasticizer_vol_fraction"]
            + output["water_vol_fraction"]
            + output["sub_vol_fraction"]
            + output["cem_vol_fraction"]
            + input_dic["aggregates_volume_fraction"]
        ).magnitude
    )

    # computation of volume contents
    output["cem_vol_content"] = output["cem_vol_fraction"] * input_dic["density_cem"]
    output["sub_vol_content"] = output["sub_vol_fraction"] * input_dic["density_sub"]
    output["water_vol_content"] = output["water_vol_fraction"] * input_dic["density_water"]
    output["aggregates_vol_content"] = input_dic["aggregates_volume_fraction"] * input_dic["density_aggregates"]

    # sanity check, the total volume adds up to 1
    assert 1 == pytest.approx(
        (
            input_dic["plasticizer_volume_content"] / input_dic["density_plasticizer"]
            + output["water_vol_content"] / input_dic["density_water"]
            + output["sub_vol_content"] / input_dic["density_sub"]
            + output["cem_vol_content"] / input_dic["density_cem"]
            + output["aggregates_vol_content"] / input_dic["density_aggregates"]
        ).magnitude
    )

    # paste density
    output["density_paste"] = (
        output["cem_vol_content"]
        + output["sub_vol_content"]
        + output["water_vol_content"]
        + input_dic["plasticizer_volume_content"]
    ) / (1 - input_dic["aggregates_volume_fraction"])

    return output


def computation_ratios(input_dic):
    """
    This is the function to compute mass or volume ratios based on volume contents
    It is basically the inverse of 'computation_volume_content'.

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    input : dic / pint units required converted to matching units
        required:
            - 'water_vol_content'
            - 'density_water'
            - 'density_cem'
            - 'cem_vol_content'
        optional (otherwise assumed to be zero)
            - 'sub_vol_content'
            - 'density_sub' : density of the substitute (slag or cemII)
            - 'aggregates_vol_content'
            - 'density_aggregates'
            - 'plasticizer_volume_content'
            - 'density_plasticizer'

    Returns
    -------
    output : dic / pint unit given
            - 'wb_mass_ratio' : water to binder ratio water/(cement + substitute (e.g. slag or cemII))
            - 'aggregates_volume_fraction'
            - 'sc_volume_fraction' : volume fraction of substitute to base cement
    """
    # optional paramters
    # plasticizer
    if "plasticizer_volume_content" not in input_dic.keys():
        input_dic["plasticizer_volume_content"] = 0 * ureg("kg/m^3")
        input_dic["density_plasticizer"] = 42 * ureg("kg/m^3")
    # substitute
    if "sub_vol_content" not in input_dic.keys():
        input_dic["sub_vol_content"] = 0 * ureg("kg/m^3")
        input_dic["density_sub"] = 42 * ureg("kg/m^3")  # dummy value
    # aggregates
    if "aggregates_vol_content" not in input_dic.keys():
        input_dic["aggregates_vol_content"] = 0 * ureg("kg/m^3")
        input_dic["density_aggregates"] = 42 * ureg("kg/m^3")  # dummy value

    # set correct units
    input_dic["density_cem"].ito("kg/m^3")
    input_dic["density_sub"].ito("kg/m^3")
    input_dic["density_water"].ito("kg/m^3")
    input_dic["density_plasticizer"].ito("kg/m^3")
    input_dic["density_aggregates"].ito("kg/m^3")

    input_dic["plasticizer_volume_content"].ito("kg/m^3")
    input_dic["cem_vol_content"].ito("kg/m^3")
    input_dic["sub_vol_content"].ito("kg/m^3")
    input_dic["water_vol_content"].ito("kg/m^3")
    input_dic["aggregates_vol_content"].ito("kg/m^3")

    output = {}

    # compute water to binder mass ratio
    ## problem: plasticizer is assued as water volume, but has different density
    ## is is no practical problem, as plastizicer is very close to the denisty of water, but numerically it is...
    ## therefore: compute plastizicer vol content as if it had density of water
    pl_as_water_content = (
        input_dic["plasticizer_volume_content"] / input_dic["density_plasticizer"] * input_dic["density_water"]
    )

    output["wb_mass_ratio"] = (input_dic["water_vol_content"] + pl_as_water_content) / (
        input_dic["cem_vol_content"] + input_dic["sub_vol_content"]
    )

    # compute substitute to base cement volume fraction
    sub_vol = input_dic["sub_vol_content"] / input_dic["density_sub"]
    cem_vol = input_dic["cem_vol_content"] / input_dic["density_cem"]

    # aggregate volume fraction
    output["aggregates_volume_fraction"] = input_dic["aggregates_vol_content"] / input_dic["density_aggregates"]

    # compute substitute to base cement ratio
    output["sc_volume_fraction"] = sub_vol / (sub_vol + cem_vol)

    return output
