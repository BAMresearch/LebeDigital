import numpy as np
from lebedigital.unit_registry import ureg

@ureg.wraps(('kg_CO2_eq'),('kg_CO2_eq/m^3', 'm', 'm', 'm'))
def computation_GWP_per_part(gwp_mix, width, height, length):
    """
    This function computes the global warming potential for a single beam

    Parameters
    ----------
    gwp_mix : float / pint unit in kg_CO2_eq/m^3
        GWP per cubic meter of used mix
    width : float / pint unit length
        width of the beam
    height : float / pint unit length
        height of the beam
    length : float / pint unit length
        length of the beam

    Returns
    -------
    beam_gwp : float / pint unit, will be in 'kg_CO2_eq'
       GWP of one beam
    """

    beam_gwp = width * height * length * gwp_mix

    return beam_gwp