import fenics_concrete

def concrete_homogenization(parameters):
    """ returns homogenized concrete parameter

    this function calls the Mori-Tanaka homogenization scheme
    A dictionary with the required parameter is given

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


    # initialize concrete paste
    concrete = fenics_concrete.ConcreteHomogenization(E_matrix=parameters['paste_E'],
                                                      nu_matrix=parameters['paste_nu'],
                                                      fc_matrix=parameters['paste_fc'],
                                                      kappa_matrix=parameters['paste_kappa'],
                                                      rho_matrix=parameters['paste_rho'],
                                                      C_matrix=parameters['paste_C'],
                                                      Q_matrix=parameters['paste_Q'])

    # adding uncoated aggregates
    concrete.add_uncoated_particle(E=parameters['aggregates_E'],
                                   nu=parameters['aggregates_nu'],
                                   volume_fraction=parameters['aggregates_vol_frac'],
                                   kappa=parameters['aggregates_kappa'],
                                   rho=parameters['aggregates_rho'],
                                   C=parameters['aggregates_C'])

    results = {'E': concrete.E_eff,
               'nu': concrete.nu_eff,
               'fc': concrete.fc_eff,
               'C': concrete.C_vol_eff,
               'rho': concrete.rho_eff,
               'kappa': concrete.kappa_eff,
               'Q': concrete.Q_vol_eff}

    return results