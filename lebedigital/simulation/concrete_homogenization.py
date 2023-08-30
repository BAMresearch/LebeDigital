import fenicsxconcrete

from lebedigital.unit_registry import ureg


def concrete_homogenization(parameters):
    """returns homogenized concrete parameter

    this function calls the Mori-Tanaka homogenization scheme
    A dictionary with the required parameter is given

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.
    Units will be automatically converted to what is required, an appropriate dimensionality must be given.

    Parameters
    ----------
    parameters : dict
        List with all required parameters
        Paste properties
        - paste_E     : Young's modulus
        - paste_nu    : Poission's ratio
        - paste_fc    : Compressive strength
        - paste_kappa : Thermal conductivity
        - paste_rho   : Density
        - paste_C     : Specific heat capacity
        - paste_Q     : Total heat release in energy per weight
        Aggregate properties
        - aggregates_E        : Young's modulus
        - aggregates_nu       : Poission's ratio
        - aggregates_vol_frac : Volume fraction
        - aggregates_kappa    : Thermal conductivity
        - aggregates_rho      : Density
        - aggregates_C        : Specific heat capacity

    Returns
    -------
    results : dict
        List with homogenized values
        - E, nu, fc, C, rho, kappa, Q
    """

    # converting to correct pint units / automatic check for pint input
    parameters["paste_E"].ito("Pa")
    parameters["paste_fc"].ito("Pa")
    parameters["paste_kappa"].ito("W/m/K")
    parameters["paste_rho"].ito("kg/m^3")
    parameters["paste_C"].ito("J/kg/K")
    parameters["paste_Q"].ito("J/kg")
    parameters["aggregates_E"].ito("Pa")
    parameters["aggregates_kappa"].ito("W/m/K")
    parameters["aggregates_rho"].ito("kg/m^3")
    parameters["aggregates_C"].ito("J/kg/K")

    # initialize concrete paste
    concrete = fenicsxconcrete.ConcreteHomogenization(
        E_matrix=parameters["paste_E"].magnitude,
        nu_matrix=parameters["paste_nu"].magnitude,
        fc_matrix=parameters["paste_fc"].magnitude,
        kappa_matrix=parameters["paste_kappa"].magnitude,
        rho_matrix=parameters["paste_rho"].magnitude,
        C_matrix=parameters["paste_C"].magnitude,
        Q_matrix=parameters["paste_Q"].magnitude,
    )

    # adding uncoated aggregates
    concrete.add_uncoated_particle(
        E=parameters["aggregates_E"].magnitude,
        nu=parameters["aggregates_nu"].magnitude,
        volume_fraction=parameters["aggregates_vol_frac"].magnitude,
        kappa=parameters["aggregates_kappa"].magnitude,
        rho=parameters["aggregates_rho"].magnitude,
        C=parameters["aggregates_C"].magnitude,
    )

    # output with corresponding units
    results = {
        "E": concrete.E_eff * ureg("Pa"),
        "nu": concrete.nu_eff * ureg("dimensionless"),
        "fc": concrete.fc_eff * ureg("Pa"),
        "C": concrete.C_vol_eff * ureg("J/m^3/K"),
        "rho": concrete.rho_eff * ureg("kg/m^3"),
        "kappa": concrete.kappa_eff * ureg("W/m/K"),
        "Q": concrete.Q_vol_eff * ureg("J/m^3"),
    }

    return results
