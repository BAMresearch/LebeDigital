from __future__ import print_function
import fenics_concrete
from lebedigital.simulation.precast_column import column_simulation
from lebedigital.simulation.concrete_homogenization import concrete_homogenization
import pandas as pd

# setting up the problem
def column_simulation_plus_homogenization(hom_params, simulation_params, simulation_settings):
        # run simulation

        homogenization_results = concrete_homogenization(hom_params)



        # concrete values
        simulation_params['density'] = homogenization_results['rho']   # in kg/m^3 density of concrete
        simulation_params['themal_cond'] = homogenization_results['kappa'] # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
        simulation_params['vol_heat_cap'] = homogenization_results['C']  # volumetric heat cap J/m3

        # simulation_params for mechanics problem
        simulation_params['E_28'] = homogenization_results['E']  # Youngs Modulus in Pa
        simulation_params['nu'] = homogenization_results['nu']  # Poissons Ratio

        simulation_params['fc_inf'] = homogenization_results['fc']  # in Pa
        simulation_params['ft_inf'] = homogenization_results['fc']/10  # in Pa APPROXIMATION!!!

        # Q_inf: computed as Q_pot (heat release in J/kg of binder) * density binder * vol_frac. of binder
        simulation_params['Q_inf'] = homogenization_results['Q'] # potential heat per volume of concrete in J/m^3

        #simulation_params['Q_inf'] = 240000000  # potential heat per volume of concrete in J/m^3

        data = column_simulation(simulation_settings['full_time'], simulation_settings['time_step'], simulation_params)

        return data

        #
        # assert data['time'].tolist() == pytest.approx([1200, 2400, 3600])
        # assert data['temperature'].tolist() == pytest.approx([41.487825, 43.581025, 48.334999])
        # assert data['yield'].tolist() == pytest.approx([129715.771538, 100205.750197, 46113.785397])



if __name__ == "__main__":
        # homogenization paramters
        # initialize dictionary
        hom_params = {}
        # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
        # kg - m - s - N  - Pa - J
        # values are kind of made up but within the expected magnitude
        # paste data
        hom_params['paste_E'] = 30e9  # Pa
        hom_params['paste_nu'] = 0.2
        hom_params['paste_C'] = 870  # J/kg Specific Heat Capacity
        hom_params['paste_kappa'] = 1.8  # W/m/K Thermal conductivity
        hom_params['paste_rho'] = 2400  # kg/m^3
        hom_params['paste_fc'] = 30e6  # Pa
        hom_params['paste_Q'] = 250  # J/kg  TODO DOUBLE CHECK VALUE!!! bzw, final Q!!!

        # aggregate data
        hom_params['aggregates_E'] = 25e9  # Pa
        hom_params['aggregates_nu'] = 0.3
        hom_params['aggregates_C'] = 840  # J/kg Specific Heat Capacity
        hom_params['aggregates_kappa'] = 0.8  # W/m/K Thermal conductivity
        hom_params['aggregates_rho'] = 2600  # kg/m^3
        hom_params['aggregates_vol_frac'] = 0.6

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
        simulation_params['B1'] = 2.916E-4  # in 1/s
        simulation_params['B2'] = 0.0024229  # -
        simulation_params['eta'] = 5.554  # something about diffusion
        simulation_params['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
        simulation_params['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
        simulation_params['T_ref'] = 25  # reference temperature in degree celsius
        
        

        # simulation time
        simulation_settings = {'full_time' : 60 * 60 * 1,
                               'time_step' : 60 * 20}

        column_simulation_plus_homogenization(hom_params,simulation_params,simulation_settings)

