from lebedigital.unit_registry import ureg


@ureg.wraps("J/kg/K", ("", "", "", "J/kg/K", "J/kg/K", "J/kg/K", "kg/m^3", "kg/m^3", "kg/m^3"))
def computation_specific_heat_capacity_paste(
    vol_frac_cement,
    vol_frac_sub,
    vol_frac_water,
    shc_cement,
    shc_sub,
    shc_water,
    density_cem,
    density_sub,
    density_water,
):
    """
    This function computes the specific heat capacity of the cement paste

    Parameters
    ----------
    vol_frac_cement :
        volume fraction of cement for 1m^3 of concrete
    vol_frac_sub :
        volume fraction of substitute (slag) for 1m^3 of concrete
    vol_frac_water :
        volume fraction of water for 1m^3 of concrete
    shc_cement :
        specific heat capacity of cement
    shc_sub :
        specific heat capacity of substitute (slag)
    shc_water :
        specific heat capacity of water
    density_cem :
        density of cement
    density_sub :
        density of substitute (slag)
    density_water :
        density of water

    Returns
    -------
    specific_heat_capacity_paste : float / pint unit
    """
    mass_cement = vol_frac_cement * density_cem
    mass_sub = vol_frac_sub * density_sub
    mass_water = vol_frac_water * density_water

    # average wrt mass
    specific_heat_capacity_paste = (mass_cement * shc_cement + mass_sub * shc_sub + mass_water * shc_water) / (
        mass_cement + mass_sub + mass_water
    )

    return specific_heat_capacity_paste
