import fenics_concrete
from lebedigital.simulation.precast_column import column_simulation
from lebedigital.simulation.precast_column_plus_homogenization import column_simulation_plus_homogenization
import numpy as np
from scipy.interpolate import interp1d

def optimization_wrapper_column_simulation_plus_homogenization(latents : list):
    """

    Args:
        latents: Indexed as 0-B1, 1 - B2, 2, eta, 3 - paste_Q, 4 - ratio_concrete

    Returns:

    """
    hom_params = {}
    # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
    # kg - m - s - N  - Pa - J
    # values are kind of made up but within the expected magnitude
    # paste data
    hom_params['paste_E'] = 25e9  # Pa
    hom_params['paste_nu'] = 0.2
    hom_params['paste_C'] = 2.4e6/2500.0  # J/kg Specific Heat Capacity
    hom_params['paste_kappa'] = 2.0  # W/m/K Thermal conductivity
    hom_params['paste_rho'] = 2500  # kg/m^3
    hom_params['paste_fc'] = 30e6  # Pa
    #hom_params['paste_Q'] = 250  # J/kg  TODO DOUBLE CHECK VALUE!!! bzw, final Q!!!
    hom_params['paste_Q'] = latents[3]
    # aggregate data - for agg_vol_frac = 0, should be irrelevant
    hom_params['aggregates_E'] = 25e9  # Pa
    hom_params['aggregates_nu'] = 0.3
    hom_params['aggregates_C'] = 840  # J/kg Specific Heat Capacity
    hom_params['aggregates_kappa'] = 0.8  # W/m/K Thermal conductivity
    hom_params['aggregates_rho'] = 2600  # kg/m^3
    #hom_params['aggregates_vol_frac'] = 0.6
    hom_params['aggregates_vol_frac'] = 1 - latents[-1] # 1- concrete_ratio
    # setup simulation parameters:
    simulation_params = fenics_concrete.Parameters()
    # model simulation_params
    # required simulation_params for alpha to E mapping
    simulation_params['alpha_t'] = 0.2
    simulation_params['alpha_0'] = 0.05
    simulation_params['a_E'] = 0.6
    # required simulation_params for alpha to tensile and compressive stiffness mapping
    simulation_params['a_fc'] = 1.5
    simulation_params['a_ft'] = 1.2
    # temperature setings:
    simulation_params['T_0'] = 40  # initial temperature of concrete
    simulation_params['T_bc1'] = 20  # constant boundary temperature
    # column geometry
    simulation_params['width'] = 0.5  # length of pillar in m
    simulation_params['height'] = 4  # width (square cross-section)
    # values for hydration
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    simulation_params['B1'] = latents[0] # 2.916E-4  # in 1/s
    simulation_params['B2'] = latents[1] # 0.0024229  # -
    simulation_params['eta'] = latents[2] # 5.554  # something about diffusion
    simulation_params['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
    simulation_params['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
    simulation_params['T_ref'] = 25  # reference temperature in degree celsius

    # simulation time
    simulation_settings = {'full_time': 60 * 60 * 5,
                           'time_step': 60 * 20}

    data = column_simulation_plus_homogenization(hom_params, simulation_params, simulation_settings)

    # --- Values specific for the optimisation problem
    # time at which the yield turns to negative, or the column can sustain its own weight.
    f = interp1d(data['yield'], data['time'], kind='cubic')
    time_critical = f([0.])  # returs the time when the yield point has reached

    # Max temp attained by the column in the simulation time. It attains a max value and then falls down
    temp_max = np.max(data['temperature'])

    return data, time_critical, temp_max

if __name__ == '__main__':
    #testing
    #scaling = np.array([1e-04,1e-03,1,1e04,1])
    #latents = np.array([2, 6.32, 3.5, 2.5,1])*scaling
    #latents = np.array([2.916E-4,0.0024229,5.554,250,0.4])

    B1 = 2E-4  # in 1/s
    B2 = 6.32e-03  # -
    eta = 3.5
    paste_Q = 4.2e5*0.2
    ratio_concrete = 1
    latents = np.array([B1, B2, eta, paste_Q,ratio_concrete])


    data,time_critical, temp_max = optimization_wrapper_column_simulation_plus_homogenization(latents)
    print(data)