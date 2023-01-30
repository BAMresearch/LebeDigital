import numpy as np
from lebedigital.unit_registry import ureg

def computation_GWP_mix(constituents):
    """
    This function computes the global warming potential for a unit of beam

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

    gwp_mix = 0 * ureg('kg_CO2_eq/m^3')
    for constituent in constituents.keys():
        # check and convert units
        constituents[constituent]['content'].ito('kg/m^3')
        constituents[constituent]['GWP'].ito('kg_CO2_eq/kg')
        gwp_mix += constituents[constituent]['content'] * constituents[constituent]['GWP']

    return gwp_mix

if __name__ == "__main__":
    # test while developing this
    constituents = {'concrete' : {'content': 5 * ureg('kg/m^3'), 'GWP': 2 * ureg('kg_CO2_eq/kg') },
                    'aggregates': {'content': 10 * ureg('kg/m^3'), 'GWP': 1 * ureg('kg_CO2_eq/kg') },
                    'slag': {'content': 5 * ureg('kg/m^3'), 'GWP': 20 * ureg('kg_CO2_eq/kg') }
                    }


    print(computation_GWP_mix(constituents))
