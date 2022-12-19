from lebedigital.unit_registry import ureg

@ureg.check('[stress]')
def approximate_tensile_strength(compressive_strength):
    """
    function to approximate tensile strength from compressive strength

    Parameters
    ----------
    compressive_strength : float
        fc for concrete

    Returns
    -------
    tensile_strength : float
        concrete compressive strength
    """
    tensile_strength = compressive_strength/10  # rough approximation

    return tensile_strength
