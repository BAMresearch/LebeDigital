import fenics_concrete
import pandas as pd
import pint_pandas
from lebedigital.unit_registry import ureg

def demonstrator_beam(time, dt, parameters, pv_output=False, pv_name='beam_simulation'):

        # check/convert units...
        # TODO this needs to be moved to fenics concrete but for now this is ok.
        parameters["density"].ito('kg/m^3')
        parameters["themal_cond"].ito('W/m/K')
        parameters['vol_heat_cap'].ito('J/m^3/K')
        parameters["alpha_t"].ito('')
        parameters["alpha_0"].ito('')
        parameters["a_E"].ito('')
        parameters["fc_inf"].ito('N/m^2')
        parameters["a_fc"].ito('')
        parameters["ft_inf"].ito('N/m^2')
        parameters["a_ft"].ito('')
        parameters["T_0"].ito('degree_Celsius')
        parameters["T_bc1"].ito('degree_Celsius')
        parameters["length"].ito('m')
        parameters["width"].ito('m')
        parameters["height"].ito('m')
        parameters["Q_inf"].ito('J/m^3')
        parameters["B1"].ito('1/s')
        parameters["B2"].ito('')
        parameters["eta"].ito('')
        parameters["E_act"].ito('J/mol')
        parameters["T_ref"].ito('degree_Celsius')
        parameters["alpha_max"].ito('')
        parameters['E_28'].ito('N/m^2')
        parameters['nu'].ito('')
        parameters['mesh_density'].ito('')
        parameters['mesh_density_min'].ito('')

        # ... now remove all units, otherwise problems with simulation :(, will change in the future...
        for key in parameters.keys():
            # check if paramter has a unit, then remove it
            if type(parameters[key]) == type(1 * ureg('')):
                parameters[key] = parameters[key].magnitude

        # same for the two values, not in the parameter list
        time.ito('s')
        time = time.magnitude
        dt.ito('s')
        dt = dt.magnitude

        # simulation parameters
        parameters['log_level'] = 'WARNING'
        parameters['dim'] = 3
        parameters['evolution_ft'] = False # just gravity
        parameters['bc_setting'] = 'no_external_load'  # just gravity

        experiment = fenics_concrete.ConcreteBeamExperiment(parameters)
        problem = fenics_concrete.ConcreteThermoMechanical(experiment, parameters, pv_name=pv_name)

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

        # Building Pandas-Pint DataFrame
        pint_df = pd.DataFrame({"time": pd.Series(problem.sensors['MaxYieldSensor'].time, dtype="pint[s]"),
                                "temperature": pd.Series(problem.sensors['MaxTemperatureSensor'].data,
                                                         dtype="pint[degree_Celsius]"),
                                "yield": pd.Series(problem.sensors['MaxYieldSensor'].data, dtype='pint[]')})

        return pint_df