from lebedigital.unit_registry import ureg
import numpy as np

@ureg.wraps('MPa', 'MPa')
def approximate_tensile_strength(compressive_strength):
    """
    function to approximate tensile strength from compressive strength
    the computation is based on DIN EN 1992-1-1

    The input is checked to be using the python pint package: https://pint.readthedocs.io/
    This requires the correct units to be attached to the input value.

    Parameters
    ----------
    compressive_strength : float / pint stress unit
        characteristic compressive stength of a concrete cylinder

    Returns
    -------
    tensile_strength : float / pint stress unit, the same unit as the input
        concrete compressive strength
    """

    if compressive_strength <= 50:
        tensile_strength = 0.3 * compressive_strength**(2/3)
    else:
        tensile_strength = 2.12 * np.log(1+((compressive_strength+8)/10))

    return tensile_strength