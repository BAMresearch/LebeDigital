from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import fenics_concrete

# setting up the problem
def column_simulation(time, dt, parameters):
    time = time * 60 * 60  # conversion to seconds
    dt = dt  * 60  # conversion to seconds

    # simulation parameters
    parameters['mesh_density'] = 5
    parameters['mesh_density_min'] = 4
    parameters['log_level'] = 'WARNING'
    parameters['dim'] = 3
    parameters['bc_setting'] = 'full'  # default boundary setting

    experiment = fenics_concrete.ConcreteColumnExperiment(parameters)
    problem = fenics_concrete.ConcreteThermoMechanical(experiment, parameters)

    problem.add_sensor(fenics_concrete.sensors.MaxYieldSensor())
    problem.add_sensor(fenics_concrete.sensors.MaxTemperatureSensor())

    # set time step
    problem.set_timestep(dt)  # for time integration scheme

    # initialize time
    t = dt  # first time step time

    print('Run simulation!')
    time_data = []
    while t <= time:  # time loop in seconds
        # solve temp-hydration-mechanics
        problem.solve(t=t)  # solving this
        #problem.pv_plot(t=t)

        # prepare next timestep
        t += dt

    print('Done!')
    # return lists with time steps, max temperature, max yield
    return problem.sensors['MaxYieldSensor'].time, problem.sensors['MaxTemperatureSensor'].data, problem.sensors['MaxYieldSensor'].data

# example of how to use this function with potential parameters
if __name__ == "__main__":

    # setup parameters:
    parameters = fenics_concrete.Parameters()

    # model parameters
    # concrete values
    parameters['density'] = 2350  # in kg/m^3 density of concrete
    parameters['density_binder'] = 1440  # in kg/m^3 density of the binder
    parameters['themal_cond'] = 2.0  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
    parameters['vol_heat_cap'] = 2.4e6  # volumetric heat cap J/(m3 K)

    # parameters for mechanics problem
    parameters['E_28'] = 15000000  # Youngs Modulus N/m2 or something... TODO: check units!
    parameters['nu'] = 0.2  # Poissons Ratio

    # required parameters for alpha to E mapping
    parameters['alpha_t'] = 0.2
    parameters['alpha_0'] = 0.05
    parameters['a_E'] = 0.6

    # required parameters for alpha to tensile and compressive stiffness mapping
    parameters['fc_inf'] = 6210000
    parameters['a_fc'] = 1.2
    parameters['ft_inf'] = 467000
    parameters['a_ft'] = 1.0

    # temperature setings:
    parameters['T_0'] = 40  # initial temperature of concrete
    parameters['T_boundary'] = 40  # constant boundary temperature

    # column geometry
    parameters['width'] = 0.5  # length of pillar in m
    parameters['height'] = 4 # width (square cross-section)

    # values for hydration
    parameters['b_ratio'] = 0.2  # volume percentage of binder
    parameters['Q_pot'] = 500e3  # potential heat per weight of binder in J/kg
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    parameters['B1'] = 2.916E-4  # in 1/s
    parameters['B2'] = 0.0024229  # -
    parameters['eta'] = 5.554  # something about diffusion
    parameters['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
    parameters['E_act'] = 5653 * 8.3145  # activation energy in Jmol^-1
    parameters['T_ref'] = 25  # reference temperature in degree celsius

    # simulation time
    full_time = 4  # simulation time in hours
    time_step = 10  # timestep in minutes

    # run simulation
    time_data, temp_data, yield_data = column_simulation(full_time, time_step, parameters)

    print('Print data:')
    for i in range(len(time_data)):
        print(time_data[i], temp_data[i], yield_data[i])

    # TODO: fix initial data for sensors!!!!
    #       fix time list for some sensors!!!
    #       check material parameters, units and values
    #       how is the cement etc. implemented...








