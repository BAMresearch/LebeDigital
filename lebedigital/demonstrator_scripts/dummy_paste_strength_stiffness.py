import numpy as np
from lebedigital.unit_registry import ureg

@ureg.check('','')
def dummy_paste_strength_stiffness(slag_ratio, phi_paste) :
    """
    This is a dummy function to make the snakemake workflow work until the real function is ready
    It changes arbitrarily chosen values depending on slag content, not based on physics or anything.

    Parameters
    ----------
    slag_ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    phi_paste: ??
        input from Atuls parameter identification, not sure currently

    Returns
    -------
    paste_youngs_modulus : float / pint stress unit, will be in 'GPa'
        approximated youngs modulus of paste
    paste_strength : float / pint stress unit, will be in 'MPa'
        approximated compressive strength of paste
    """

    paste_youngs_modulus_min = 15 * ureg('GPa')
    paste_youngs_modulus_max = 25 * ureg('GPa')
    paste_youngs_modulus = paste_youngs_modulus_min + (paste_youngs_modulus_max-paste_youngs_modulus_min)*slag_ratio

    paste_strength_min = 8 * ureg('MPa')
    paste_strength_max = 13 * ureg('MPa')
    paste_strength = paste_strength_min + (paste_strength_max-paste_strength_min)*slag_ratio

    return paste_youngs_modulus, paste_strength

if __name__ == "__main__":
    # test while developing this
    E, fc = dummy_paste_strength_stiffness(0,10)
    print(E, fc)
