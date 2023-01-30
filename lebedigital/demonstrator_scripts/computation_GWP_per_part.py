import numpy as np
from lebedigital.unit_registry import ureg

@ureg.wraps(('kg_CO2_eq'),('kg_CO2_eq/m^3', 'm', 'm', 'm'))
def computation_GWP_per_part(gwp, width, height, length):
    """
    This is a dummy function to make the snakemake workflow work until the real function is

    Parameters
    ----------
    slag_ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    phi_hydration: ??
        input from Atuls parameter identification

    Returns
    -------
    B1 : float / pint unit will be in '1/s'
        hydration parameter
    B2 : float / pint unitless unit
        hydration parameter
    eta : float / pint unitless unit
        hydration parameter
    E_act : float / pint unit, will be in 'J/mol'
        activation energy - hydration parameter
    T_ref : float / pint unit, will be in 'degree Celsius'
        reference temperature - hydration parameter
    Q_pot : float / pint unit, will be in 'J/kg'
        maximum potential hydration parameter
    """

    return width * height * length *gwp

if __name__ == "__main__":
    # test while developing this
    gwp = 10 * ureg('kg_CO2_eq/m^3')
    width = 0.2 * ureg('m')
    height = 0.5 * ureg('m')
    length = 10000 * ureg('mm')


    print(computation_GWP_per_part(gwp,width,height,length))
