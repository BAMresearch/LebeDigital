import pytest
import fenics_concrete
from lebedigital.simulation.precast_column import column_simulation

def test_column_simulation():
    # setup parameters:
    parameters = fenics_concrete.Parameters()

    # model parameters
    # concrete values
    parameters['density'] = 2350  # in kg/m^3 density of concrete
    parameters['themal_cond'] = 2.0  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
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
    parameters['ft_inf'] = 4e6  # in Pa
    parameters['a_ft'] = 1.2

    # temperature setings:
    parameters['T_0'] = 40  # initial temperature of concrete
    parameters['T_bc1'] = 20  # constant boundary temperature

    # column geometry
    parameters['width'] = 0.5  # length of pillar in m
    parameters['height'] = 4  # width (square cross-section)

    # values for hydration
    # Q_inf: computed as Q_pot (heat release in J/kg of binder) * density binder * vol_frac. of binder
    parameters['Q_inf'] = 240000000  # potential heat per volume of concrete in J/m^3
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    parameters['B1'] = 2.916E-4  # in 1/s
    parameters['B2'] = 0.0024229  # -
    parameters['eta'] = 5.554  # something about diffusion
    parameters['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
    parameters['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
    parameters['T_ref'] = 25  # reference temperature in degree celsius

    # simulation time
    full_time = 60*60*1  # simulation time in hours
    time_step = 60*20  # timestep in minutes

    # run simulation
    data = column_simulation(full_time, time_step, parameters)

    assert data['time'].tolist() == pytest.approx([1200, 2400, 3600])
    assert data['temperature'].tolist() == pytest.approx([41.487825, 43.581025, 48.334999])
    assert data['yield'].tolist() == pytest.approx([129715.771538, 100205.750197, 46113.785397])