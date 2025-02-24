import pytest
import fenics_concrete
from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.simulation.precast_column_plus_homogenization import column_simulation_plus_homogenization
from lebedigital.simulation.precast_column import column_simulation


def test_column_plus_homogenization():
    # making sure that


    #### only column
    # setup parameters:
    parameters = fenics_concrete.Parameters()

    # model parameters
    # concrete values
    parameters['density'] = 2500.0  # in kg/m^3 density of concrete
    parameters['themal_cond'] = 2  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
    parameters['vol_heat_cap'] = 2.4e6  # volumetric heat cap J/m3

    # parameters for mechanics problem
    parameters['E_28'] = 25e9  # Youngs Modulus in Pa
    parameters['nu'] = 0.2  # Poissons Ratio

    # required parameters for alpha to E mapping
    parameters['alpha_t'] = 0.2
    parameters['alpha_0'] = 0.05
    parameters['a_E'] = 0.6

    # required parameters for alpha to tensile and compressive stiffness mapping
    parameters['fc_inf'] = 30e6  # in Pa
    parameters['a_fc'] = 1.5
    parameters['ft_inf'] = parameters['fc_inf']/10  # in Pa
    parameters['a_ft'] = 1.2

    # temperature setings:
    parameters['T_0'] = 40  # initial temperature of concrete
    parameters['T_bc1'] = 20  # constant boundary temperature

    # column geometry
    parameters['width'] = 0.5  # length of pillar in m
    parameters['height'] = 4  # width (square cross-section)

    # values for hydration
    # Q_inf: computed as Q_pot (heat release in J/kg of binder) * density binder * vol_frac. of binder
    # Choose something, take 2500 kg/m³ as density and maybe something between 0.3 and 0.5 as volume fraction
    # (needs to be > 0 and <= 1). The vol fraction is basically a possibility to increase or reduce your heat output.
    # So if for a vol_frac of 0.5 your temperature exceeds your limit, you could reduce your amount of cement
    # (the thing that generates the heat).

    parameters['Q_inf'] = 4.2e5*0.2*parameters['density']  # potential heat per volume of concrete in J/m^3 ## 240000000
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    parameters['B1'] = 2E-4  # in 1/s ## 1.5 * 2.916E-4  # in 1/s
    parameters['B2'] = 6.32e-03  # - ## 0.0024229
    parameters['eta'] = 3.5  # something about diffusion  ## 5.554
    parameters['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
    parameters['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
    parameters['T_ref'] = 25  # reference temperature in degree celsius

    # simulation time
    full_time = 60 * 60 * 5  # simulation time in hours
    time_step = 60 * 20  # timestep in minutes

    # run simulation
    c_data = column_simulation(full_time, time_step, parameters)

############# c + h
    # homogenization paramters
    # initialize dictionary
    hom_params = {}
    # using consistent units  https://www.dynasupport.com/howtos/general/consistent-units
    # kg - m - s - N  - Pa - J
    # values are kind of made up but within the expected magnitude
    # paste data
    hom_params['paste_E'] = parameters['E_28']  # Pa
    hom_params['paste_nu'] = parameters['nu']
    hom_params['paste_C'] = parameters['vol_heat_cap']/parameters['density']  # J/kg Specific Heat Capacity        parameters['vol_heat_cap']????
    hom_params['paste_kappa'] = parameters['themal_cond']  # W/m/K Thermal conductivity
    hom_params['paste_rho'] = parameters['density']  # kg/m^3
    hom_params['paste_fc'] = parameters['fc_inf']  # Pa
    hom_params['paste_Q'] = parameters['Q_inf']/parameters['density']  # J/kg


    # aggregate data # should all be irrelevant with vol = 0
    hom_params['aggregates_E'] = 25e9  # Pa
    hom_params['aggregates_nu'] = 0.3
    hom_params['aggregates_C'] = 840  # J/kg Specific Heat Capacity
    hom_params['aggregates_kappa'] = 0.8  # W/m/K Thermal conductivity
    hom_params['aggregates_rho'] = 2600  # kg/m^3
    hom_params['aggregates_vol_frac'] = 0.0

    # setup simulation parameters:
    simulation_params = fenics_concrete.Parameters()

    # model simulation_params

    # required simulation_params for alpha to E mapping
    simulation_params['alpha_t'] = parameters['alpha_t']
    simulation_params['alpha_0'] = parameters['alpha_0']
    simulation_params['a_E'] = parameters['a_E']


    # required simulation_params for alpha to tensile and compressive stiffness mapping
    simulation_params['a_fc'] = parameters['a_fc']
    simulation_params['a_ft'] = parameters['a_ft']

    # temperature setings:
    simulation_params['T_0'] = parameters['T_0']  # initial temperature of concrete
    simulation_params['T_bc1'] = parameters['T_bc1']  # constant boundary temperature

    # column geometry
    simulation_params['width'] = parameters['width']  # length of pillar in m
    simulation_params['height'] = parameters['height']  # width (square cross-section)

    # values for hydration
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    simulation_params['B1'] = parameters['B1']  # in 1/s
    simulation_params['B2'] = parameters['B2']  # -
    simulation_params['eta'] = parameters['eta']  # something about diffusion
    simulation_params['alpha_max'] = parameters['alpha_max']  # also possible to approximate based on equation with w/c
    simulation_params['E_act'] = parameters['E_act']  # activation energy in Jmol^-1
    simulation_params['T_ref'] = parameters['T_ref']  # reference temperature in degree celsius

    # simulation time
    simulation_settings = {'full_time': full_time,
                           'time_step': time_step}

    c_plus_h_data = column_simulation_plus_homogenization(hom_params, simulation_params, simulation_settings)

    assert c_data['time'].tolist() == c_plus_h_data['time'].tolist()  # just a sanity test
    assert c_data['temperature'].tolist() == pytest.approx(c_plus_h_data['temperature'].tolist())
    assert c_data['yield'].tolist() == pytest.approx(c_plus_h_data['yield'].tolist())

    print(c_data['temperature'])




if __name__ == "__main__":
    test_column_plus_homogenization()
