from __future__ import print_function
import fenics_concrete
import pandas as pd

# setting up the problem
def column_simulation(time, dt, parameters, pv_output=False):

    # simulation parameters
    parameters['mesh_density'] = 5
    parameters['mesh_density_min'] = 4
    parameters['log_level'] = 'WARNING'
    parameters['dim'] = 3
    parameters['bc_setting'] = 'full'  # default boundary setting

    sortednames = sorted(parameters.keys(), key=lambda x: x.lower())
    print('-----------------')
    for i in sortednames:
        values = parameters[i]
        print(i, values)


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

    # Building DataFrame
    df = pd.DataFrame(list(zip(problem.sensors['MaxYieldSensor'].time,
                               problem.sensors['MaxTemperatureSensor'].data,
                               problem.sensors['MaxYieldSensor'].data)),
                      columns=['time', 'temperature', 'yield'])

    # return lists with time steps, max temperature, max yield
    return df
