import pytest
import fenics_concrete
from lebedigital.simulation.precast_column import column_simulation
import numpy as np
from scipy.interpolate import interp1d

def Column_simulation(latents :list):
    """

    Args:
        latents (): idx 0 - B_1, idx 1 - B2, idx 2 - \eta, idx 3 - Q_pot
        Note whatever scaling needs to be done, it needs to be externally done
    Returns:

    """
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

    # temperature settings:
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
    densityBinder = 2500 #kg/m3
    vol_frac_binder = 0.2
    Q_inf = latents[-1]*densityBinder*vol_frac_binder
    #parameters['Q_inf'] = 240000000  # potential heat per volume of concrete in J/m^3
    parameters['Q_inf'] = Q_inf
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    parameters['B1'] = latents[0]  # in 1/s
    parameters['B2'] = latents[1] # -
    parameters['eta'] = latents[2]  # something about diffusion
    parameters['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
    parameters['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
    parameters['T_ref'] = 25  # reference temperature in degree celsius

    # simulation time
    full_time = 60*60*5  # simulation time in hours
    time_step = 60*20  # timestep in minutes

    # run simulation
    data = column_simulation(full_time, time_step, parameters)

    # --- Values specific for the optimisation problem
    # time at which the yield turns to negative, or the column can sustain its own weight.
    f = interp1d(data['yield'],data['time'],kind='cubic')
    time_critical = f([0.]) # returs the time when the yield point has reached

    # Max temp attained by the column in the simulation time. It attains a max value and then falls down
    temp_max = np.max(data['temperature'])

    return data, time_critical, temp_max

# testing
#scaling = np.array([1e-04,1e-03,1,1e05])
#latents = np.array([2, 6.32, 3.5, 4.2])*scaling

#data,time_critical, temp_max = Column_simulation(latents)
