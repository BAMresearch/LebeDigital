import numpy as np

from lebedigital.unit_registry import ureg


@ureg.wraps("kg_CO2_eq", ("kg_CO2_eq/m^3", "kg_CO2_eq/m^3", "m", "m", "m", "", "m"))
def computation_GWP_per_part(gwp_mix, gwp_steel, width, height, length, n_steel, diameter_steel):
    """
    This function computes the global warming potential for a single beam

    Parameters
    ----------
    gwp_mix : float / pint unit in kg_CO2_eq/m^3
        GWP per cubic meter of used mix
    gwp_steel : float / pint unit in kg_CO2_eq/m^3
        GWP per cubic meter of used reinforcement steel
    width : float / pint unit length
        width of the beam
    height : float / pint unit length
        height of the beam
    length : float / pint unit length
        length of the beam
    n_steel :
        number of used steel rebars
    diameter_steel :
        diameter of used steel rebars

    Returns
    -------
    beam_gwp : float / pint unit, will be in 'kg_CO2_eq'
       GWP of one beam
    """

    # concrete
    # removing the volume of the reinforcement form the concrete calculation is ignored
    beam_gwp = width * height * length * gwp_mix
    # steel
    # difference between rebar length and beam length is ignored
    beam_gwp += n_steel * (np.pi * (diameter_steel / 2) ** 2) * length * gwp_steel

    return beam_gwp
