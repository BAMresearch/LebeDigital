import pytest
from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg
import numpy as np
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

def create_homogenization_figure(parameter_with_units : dict):
    # initialize dictionary
    parameters = {}

    # paste data
    parameters['paste_E'] = parameter_with_units['paste_E'] * ureg(parameter_with_units['paste_E_unit'])
    parameters['paste_nu'] = parameter_with_units['paste_nu'] * ureg(parameter_with_units['paste_nu_unit'])
    parameters['paste_C'] = parameter_with_units['paste_C'] * ureg(parameter_with_units['paste_C_unit'])
    parameters['paste_kappa'] = parameter_with_units['paste_kappa'] * ureg(parameter_with_units['paste_kappa_unit'])
    parameters['paste_rho'] = parameter_with_units['paste_rho'] * ureg(parameter_with_units['paste_rho_unit'])
    parameters['paste_fc'] = parameter_with_units['paste_fc'] * ureg(parameter_with_units['paste_fc_unit'])
    parameters['paste_Q'] = parameter_with_units['paste_Q'] * ureg(parameter_with_units['paste_Q_unit'])

    # aggregate data
    parameters['aggregates_E'] = parameter_with_units['aggregates_E'] * ureg(parameter_with_units['aggregates_E_unit'])
    parameters['aggregates_nu'] = parameter_with_units['aggregates_nu'] * ureg(parameter_with_units['aggregates_nu_unit'])
    parameters['aggregates_C'] = parameter_with_units['aggregates_C'] * ureg(parameter_with_units['aggregates_C_unit'])
    parameters['aggregates_kappa'] = parameter_with_units['aggregates_kappa'] * ureg(parameter_with_units['aggregates_kappa_unit'])
    parameters['aggregates_rho'] = parameter_with_units['aggregates_rho'] * ureg(parameter_with_units['aggregates_rho_unit'])

    n_points = 2
    step_size = 1/n_points
    aggregates_vol_frac_list = np.arange(0,1+step_size,step_size)
    assert aggregates_vol_frac_list.max() <= 1.0

    E_list = []
    nu_list = []
    fc_list = []
    C_list = []
    rho_list = []
    kappa_list = []
    Q_list = []

    for value in aggregates_vol_frac_list:

        parameters['aggregates_vol_frac'] = value * ureg('dimensionless')
        results = concrete_homogenization(parameters)
        E_list.append(results['E'].magnitude)
        nu_list.append(results['nu'].magnitude)
        fc_list.append(results['fc'].magnitude)
        C_list.append(results['C'].magnitude)
        rho_list.append(results['rho'].magnitude)
        kappa_list.append(results['kappa'].magnitude)
        Q_list.append(results['Q'].magnitude)

    print(C_list)
    # TODO: make a nice plot, save as image (or list and use tiks...)



if __name__ == "__main__":
    homogenization_example_parameters = {
        'paste_E': 30e9,
        'paste_E_unit': 'Pa',
        'paste_nu': 0.2,
        'paste_nu_unit': 'dimensionless',
        'paste_C': 870,
        'paste_C_unit': 'J/kg/K',
        'paste_kappa': 1.8,
        'paste_kappa_unit': 'W/m/K',
        'paste_rho': 2400,
        'paste_rho_unit': 'kg/m^3',
        'paste_fc': 30e6,
        'paste_fc_unit': 'Pa',
        'paste_Q': 250000,
        'paste_Q_unit': 'J/kg',
        'aggregates_E': 25e9,
        'aggregates_E_unit': 'Pa',
        'aggregates_nu': 0.3,
        'aggregates_nu_unit': 'dimensionless',
        'aggregates_C': 840,
        'aggregates_C_unit': 'J/kg/K',
        'aggregates_kappa': 0.8,
        'aggregates_kappa_unit': 'W/m/K',
        'aggregates_rho': 2600,
        'aggregates_rho_unit': 'kg/m^3'
    }

    create_homogenization_figure(homogenization_example_parameters)
