from lebedigital.unit_registry import ureg


@ureg.wraps("kg_CO2_eq/m^3", ("kg_CO2_eq/kg", "kg/m^3"))
def computation_GWP_steel_per_volume(gwp_steel, density_steel):
    """
    This function computes the global warming potential for steel per volume

    Parameters
    ----------
    gwp_steel : float / pint unit in kg_CO2_eq/kg
        GWP per kilogram of steel
    density_steel : float / pint unit in kg/m^3
        density of steel

    Returns
    -------
    gwp_steel_per_volume : float / pint unit, will be in 'kg_CO2_eq/m^3'
    """

    gwp_steel_per_volume = gwp_steel * density_steel

    return gwp_steel_per_volume
