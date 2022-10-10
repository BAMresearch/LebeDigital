from __future__ import print_function
import fenics_concrete

# setting up the problem
def column_simulation(time, dt, parameters, pv_output=False):
    # conversion to seconds
    time = time * 60 * 60  # hours to seconds
    dt = dt * 60  # minutes to seconds

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

        if pv_output:
            problem.pv_plot(t=t)

        # prepare next timestep
        t += dt

    print('Done!')
    # return lists with time steps, max temperature, max yield
    return problem.sensors['MaxYieldSensor'].time, problem.sensors['MaxTemperatureSensor'].data, problem.sensors['MaxYieldSensor'].data