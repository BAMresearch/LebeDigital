import numpy as np
from lebedigital.unit_registry import ureg

def computation_GWP_mix(constituents):
    """
    This function computes the global warming potential for a cubic meter of concrete

    Parameters
    ----------
    constituents : dict
        the dictonary expects e.g.:
        {'cement' : {'content': 5 * ureg('kg/m^3'), 'GWP': 2 * ureg('kg_CO2_eq/kg')},
         'aggregates': {'content': 10 * ureg('kg/m^3'), 'GWP': 1 * ureg('kg_CO2_eq/kg')}}

    Returns
    -------
    gwp_mix : float / pint unit, will be in 'kg_CO2_eq/m^3'
       GWP of one cubic meter of concrete
    """

    # initialize mix gwp, with unit
    gwp_mix = 0 * ureg('kg_CO2_eq/m^3')

    for constituent in constituents.keys():
        # check and convert units
        constituents[constituent]['content'].ito('kg/m^3')
        constituents[constituent]['GWP'].ito('kg_CO2_eq/kg')
        gwp_mix += constituents[constituent]['content'] * constituents[constituent]['GWP']

    return gwp_mix