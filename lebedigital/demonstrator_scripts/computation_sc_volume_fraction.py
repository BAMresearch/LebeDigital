import pytest

from lebedigital.unit_registry import ureg


def computation_sc_volume_fraction(input_dic):
    """
    Compute the volume fraction of the substitute (slag/mix2),
    based on the mass fraction of the substitute (slag/mix2) and the densities of the substitute and cement

    Args:
        input_dic: dictionary with the following keys:
            sc_mass_fraction: mass fraction of the substitute (slag/mix2)
            density_sub: density of the substitute (slag/mix2)
            density_cem: density of the cement

    Returns:
        sc_volume_fraction: volume fraction of the substitute (slag/mix2)
    """

    if input_dic["sc_mass_fraction"] == 0:
        return 0
    else:
        mass_sub = input_dic["sc_mass_fraction"] * ureg("kg")
        mass_cem = 1 * ureg("kg") - mass_sub
        vol_sub = mass_sub / input_dic["density_sub"]
        vol_cem = mass_cem / input_dic["density_cem"]
        sc_volume_fraction = vol_sub / (vol_sub + vol_cem)

        return sc_volume_fraction
