from lebedigital.unit_registry import ureg


@ureg.wraps("J/kg/K", ("", "", "", "J/kg/K", "J/kg/K", "J/kg/K"))
def computation_specific_heat_capacity_paste(
    vol_frac_cement, vol_frac_sub, vol_frac_water, shc_cement, shc_sub, shc_water
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

    Returns
    -------
    specific_heat_capacity_paste : float / pint unit, will be in 'kg_CO2_eq/m^3'
    """

    specific_heat_capacity_paste = (
        vol_frac_cement * shc_cement + vol_frac_sub * shc_sub + vol_frac_water * shc_water
    ) / (vol_frac_cement + vol_frac_sub + vol_frac_water)

    return specific_heat_capacity_paste
