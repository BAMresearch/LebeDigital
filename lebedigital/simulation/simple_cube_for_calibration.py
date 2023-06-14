import fenics_concrete
import numpy as np
import pandas as pd
import pint_pandas

from lebedigital.unit_registry import ureg


def setup_simple_cube(time, dt, parameters, pv_output=False):
    # check/convert units...
    parameters["edge_length"].ito("m")  # length of pillar in m
    parameters["height"] = parameters["edge_length"]
    parameters["width"] = parameters["edge_length"]

    parameters["density"].ito("kg/m^3")
    parameters["themal_cond"].ito("W/m/K")
    parameters["vol_heat_cap"].ito("J/m^3/K")
    parameters["alpha_t"].ito("")
    parameters["alpha_0"].ito("")
    parameters["a_E"].ito("")
    parameters["fc"].ito("N/m^2")
    parameters["a_fc"].ito("")
    parameters["ft"].ito("N/m^2")
    parameters["a_ft"].ito("")
    parameters["ambient_temperature"].ito("degree_Celsius")
    parameters["T_0"] = parameters["ambient_temperature"]
    parameters["T_bc1"] = parameters["ambient_temperature"]
    parameters["width"].ito("m")
    parameters["height"].ito("m")
    parameters["Q_inf"].ito("J/m^3")
    parameters["B1"].ito("1/s")
    parameters["B2"].ito("")
    parameters["eta"].ito("")
    parameters["E_act"].ito("J/mol")
    parameters["T_ref"].ito("degree_Celsius")
    parameters["alpha_max"].ito("")
    parameters["E"].ito("N/m^2")
    parameters["nu"].ito("")

    # ... now remove all units, otherwise problems with simulation :(, will change in the future...
    for key in parameters.keys():
        # check if paramter has a unit, then remove it
        if type(parameters[key]) == type(1 * ureg("")):
            parameters[key] = parameters[key].magnitude

    # same for the two values, not in the parameter list
    time.ito("s")
    time = time.magnitude
    dt.ito("s")
    dt = dt.magnitude

    # simulation parameters
    parameters["mesh_density"] = 1
    parameters["mesh_density_min"] = 1
    parameters["log_level"] = "WARNING"
    parameters["dim"] = 3
    parameters["bc_setting"] = "full"  # default boundary setting

    experiment = fenics_concrete.ConcreteColumnExperiment(parameters)
    problem = fenics_concrete.ConcreteThermoMechanical(experiment, parameters)

    sensor_location = (parameters["edge_length"] / 2, parameters["edge_length"] / 2, parameters["edge_length"] / 2)
    E_sensor = fenics_concrete.sensors.YoungsModulusSensor(sensor_location)
    fc_sensor = fenics_concrete.sensors.CompressiveStrengthSensor(sensor_location)
    doh_sensor = fenics_concrete.sensors.DOHSensor(sensor_location)

    problem.add_sensor(E_sensor)
    problem.add_sensor(fc_sensor)
    problem.add_sensor(doh_sensor)

    # set time step
    problem.set_timestep(dt)  # for time integration scheme

    return problem


def get_doh_at_28day(p, ambient_temperature=ureg.Quantity(20, ureg.degC)):
    Q_ = ureg.Quantity
    parameters = fenics_concrete.Parameters()
    parameters = parameters + p

    # model parameters
    parameters["ft"] = 4e6 * ureg("N/m^2")  # in Pa  # dummy value
    parameters["a_ft"] = 1.2 * ureg("")  # dummy value

    # temperature setings:
    parameters["ambient_temperature"] = ambient_temperature

    # column geometry
    parameters["edge_length"] = 0.01 * ureg("m")

    # simulation time
    full_time = 28 * ureg("days")  # simulation time
    time_step = 1 * ureg("hour")  # timestep

    fem_problem = setup_simple_cube(full_time, time_step, parameters)

    # same for the two values, not in the parameter list
    full_time.ito("s")
    time = full_time.magnitude
    time_step.ito("s")
    dt = time_step.magnitude

    # initialize time
    t = dt  # first time step time

    print("Run simulation!")
    time_data = []
    while t <= time:  # time loop in seconds
        # solve temp-hydration-mechanics
        fem_problem.solve(t=t)  # solving this

        # prepare next timestep
        t += dt

    print("Done!")

    return fem_problem.sensors["DOHSensor"].data[-1]  # doh_at_28 days


def get_E_and_fc_over_time(p, time_list, time_step=1 * ureg("hour"), ambient_temperature=ureg.Quantity(20, ureg.degC)):
    Q_ = ureg.Quantity
    parameters = fenics_concrete.Parameters()

    parameters = parameters + p

    # model parameters
    parameters["ft"] = 4e6 * ureg("N/m^2")  # in Pa  # dummy value
    parameters["a_ft"] = 1.2 * ureg("")  # dummy value

    # temperature setings:
    parameters["ambient_temperature"] = ambient_temperature  # ambient temperature

    # column geometry
    parameters["edge_length"] = 0.01 * ureg("m")  # length of pillar in m

    # convert all entries in time list to seconds
    for i in range(len(time_list)):
        time_list[i] = time_list[i].to_base_units()

    # sorting time list by size, converting to base unit
    time_list = sorted(time_list, key=lambda x: x.magnitude)

    full_time = time_list[-1]  # simulation time (max time value)
    min_timestep = time_list[0]

    if min_timestep < time_step:
        time_step = min_timestep

    fem_problem = setup_simple_cube(full_time, time_step, parameters)

    dt = time_step.magnitude
    time = full_time.magnitude + dt

    # initialize time
    t = dt  # first time step time

    while t <= time:  # time loop in seconds
        # solve temp-hydration-mechanics
        fem_problem.solve(t=t)  # solving this

        # prepare next timestep
        t += dt

    # Building Pandas DataFrame
    df = pd.DataFrame(
        {
            "time": pd.Series(fem_problem.sensors["YoungsModulusSensor"].time, dtype="float"),
            "E": pd.Series(fem_problem.sensors["YoungsModulusSensor"].data, dtype="float"),
            # datatype float
            "fc": pd.Series(fem_problem.sensors["CompressiveStrengthSensor"].data, dtype="float"),
        }
    )

    # check if all time values are in the list
    for t in time_list:
        if t.magnitude not in df["time"].unique():
            # add nan values if not in df

            # adding new row with zero yield
            new_line = pd.DataFrame({("time"): [t.magnitude], ("E"): [np.nan], ("fc"): [np.nan]})
            df = pd.concat([df, new_line], ignore_index=True)

    # sort table
    df = df.sort_values(by=[("time")])

    # interpolate missing values
    df = df.interpolate(method="linear", limit_direction="forward")

    # add units to df
    df_units = pd.DataFrame(
        {
            "time": pd.Series(df["time"], dtype="pint[s]"),
            "E": pd.Series(df["E"], dtype="pint[N/m^2]"),
            # datatype float
            "fc": pd.Series(df["fc"], dtype="pint[N/m^2]"),
        }
    )

    return df_units
