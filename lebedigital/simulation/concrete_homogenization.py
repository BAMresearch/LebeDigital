from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import fenics_concrete

def concrete_homogenization(parameters):
    # todo describe function

    # initialize concrete paste
    concrete = fenics_concrete.ConcreteHomogenization(E_matrix=parameters['paste_E'],
                                                      nu_matrix=parameters['paste_nu'],
                                                      fc_matrix=parameters['paste_fc'],
                                                      kappa_matrix=parameters['paste_kappa'],
                                                      rho_matrix=parameters['paste_rho'],
                                                      C_matrix=parameters['paste_C'])


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
               'C': concrete.C_eff,
               'rho': concrete.rho_eff}

    return results




# example of how to use this function with potential parameters
if __name__ == "__main__":
    print('Moin')

    # initialize dictionary
    parameters = {}

    # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
    # kg - m - s - N  - Pa - J

    # paste data
    parameters['paste_E'] = 3e11  # Pa
    parameters['paste_nu'] = 0.2
    parameters['paste_C'] = 870  # J/kg Specific Heat / Heat Capacity J/kg/K??
    parameters['paste_kappa'] = 1.8  # W/m/K Thermal conductivity
    parameters['paste_rho'] = 2400  # kg/m^3
    parameters['paste_fc'] = 3e4  # Pa

    # aggregate data
    parameters['aggregates_E'] = 2.5e11  # Pa
    parameters['aggregates_nu'] = 0.3
    parameters['aggregates_C'] = 840  # J/kg Specific Heat
    parameters['aggregates_kappa'] = 0.8  # W/m/K Thermal conductivity
    parameters['aggregates_rho'] = 2600  # kg/m^3
    parameters['aggregates_vol_frac'] = 0.6

    results = concrete_homogenization(parameters)
    print(results)
