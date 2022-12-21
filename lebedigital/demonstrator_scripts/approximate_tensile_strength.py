from lebedigital.unit_registry import ureg

@ureg.check('[stress]')
def approximate_tensile_strength(compressive_strength):
    """
    function to approximate tensile strength from compressive strength

    The input is checked to be using the python pint package: https://pint.readthedocs.io/
    This requires the correct units to be attached to the input value.

    Parameters
    ----------
    compressive_strength : float / pint stress unit
        fc for concrete

    Returns
    -------
    tensile_strength : float / pint stress unit, the same unit as the input
        concrete compressive strength
    """
    tensile_strength = compressive_strength/10  # rough approximation

    return tensile_strength