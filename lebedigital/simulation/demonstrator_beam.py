import fenics_concrete
import pandas as pd
import pint_pandas

from lebedigital.unit_registry import ureg


def demonstrator_beam(time, dt, parameters, pv_output=False, pv_name="beam_simulation"):
    # check/convert units...
    # TODO this needs to be moved to fenics concrete but for now this is ok.
    parameters["density"].ito("kg/m^3")
    parameters["themal_cond"].ito("W/m/K")
    parameters["vol_heat_cap"].ito("J/m^3/K")
    parameters["alpha_t"].ito("")
    parameters["alpha_0"].ito("")
    parameters["a_E"].ito("")
    parameters["fc_inf"].ito("N/m^2")
    parameters["a_fc"].ito("")
    parameters["ft_inf"].ito("N/m^2")
    parameters["a_ft"].ito("")
    parameters["T_0"].ito("degree_Celsius")
    parameters["T_bc1"].ito("degree_Celsius")
    parameters["length"].ito("m")
    parameters["width"].ito("m")
    parameters["height"].ito("m")
    parameters["Q_inf"].ito("J/m^3")
    parameters["B1"].ito("1/s")
    parameters["B2"].ito("")
    parameters["eta"].ito("")
    parameters["E_act"].ito("J/mol")
    parameters["T_ref"].ito("degree_Celsius")
    parameters["alpha_max"].ito("")
    parameters["E_28"].ito("N/m^2")
    parameters["nu"].ito("")
    parameters["mesh_density"].ito("")
    parameters["mesh_density_min"].ito("")

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
    parameters["log_level"] = "WARNING"
    parameters["dim"] = 3
    parameters["evolution_ft"] = False  # just gravity
    parameters["bc_setting"] = "no_external_load"  # just gravity

    experiment = fenics_concrete.ConcreteBeamExperiment(parameters)
    problem = fenics_concrete.ConcreteThermoMechanical(experiment, parameters, pv_name=pv_name)

    problem.add_sensor(fenics_concrete.sensors.MaxYieldSensor())
    problem.add_sensor(fenics_concrete.sensors.MaxTemperatureSensor())

    # set time step
    problem.set_timestep(dt)  # for time integration scheme

    # initialize time
    t = dt  # first time step time

    print("Run simulation!")
    time_data = []
    while t <= time:  # time loop in seconds
        # solve temp-hydration-mechanics
        problem.solve(t=t)  # solving this

        if pv_output:
            problem.pv_plot(t=t)

        # prepare next timestep
        t += dt

    print("Done!")

    # Building Pandas-Pint DataFrame
    pint_df = pd.DataFrame(
        {
            "time": pd.Series(problem.sensors["MaxYieldSensor"].time, dtype="pint[s]"),
            "temperature": pd.Series(problem.sensors["MaxTemperatureSensor"].data, dtype="pint[degree_Celsius]"),
            "yield": pd.Series(problem.sensors["MaxYieldSensor"].data, dtype="pint[]"),
        }
    )

    return pint_df


if __name__ == "__main__":
    parameters = fenics_concrete.Parameters()

    # model parameters

    Q_ = ureg.Quantity
    # concrete values
    parameters["density"] = 2350 * ureg("kg/m^3")  # in kg/m^3 density of concrete
    parameters["themal_cond"] = 2.0 * ureg("W/m/K")  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
    parameters["vol_heat_cap"] = 2.4e6 * ureg("J/m^3/K")  # volumetric heat cap J/m3

    # parameters for mechanics problem
    parameters["E_28"] = 25e9 * ureg("N/m^2")  # Youngs Modulus in Pa
    parameters["nu"] = 0.2 * ureg("")  # Poissons Ratio

    # required parameters for alpha to E mapping
    parameters["alpha_t"] = 0.2 * ureg("")
    parameters["alpha_0"] = 0.05 * ureg("")
    parameters["a_E"] = 0.6 * ureg("")

    # required parameters for alpha to tensile and compressive stiffness mapping
    parameters["fc_inf"] = 30e6 * ureg("N/m^2")  # in Pa
    parameters["a_fc"] = 1.5 * ureg("")
    parameters["ft_inf"] = 760 * ureg("MPa")  # in Pa
    parameters["a_ft"] = 1.2 * ureg("")

    # temperature setings:
    parameters["T_0"] = Q_(40, ureg.degC)  # initial temperature of concrete
    parameters["T_bc1"] = Q_(20, ureg.degC)  # constant boundary temperature

    # column geometry
    parameters["width"] = 0.3 * ureg("m")  # length of pillar in m
    parameters["height"] = 0.1 * ureg("m")  # width (square cross-section)
    parameters["length"] = 10 * ureg("m")  # width (square cross-section)

    # values for hydration
    # Q_inf: computed as Q_pot (heat release in J/kg of binder) * density binder * vol_frac. of binder
    parameters["Q_inf"] = 240000000 * ureg("J/m^3")  # potential heat per volume of concrete in J/m^3
    # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
    parameters["B1"] = 2.916e-4 * ureg("1/s")  # in 1/s
    parameters["B2"] = 0.0024229 * ureg("")  # -
    parameters["eta"] = 5.554 * ureg("")  # something about diffusion
    parameters["alpha_max"] = 0.875 * ureg("")  # also possible to approximate based on equation with w/c
    parameters["E_act"] = 5653 * 8.3145 * ureg("J/mol")  # activation energy in Jmol^-1
    parameters["T_ref"] = Q_(25, ureg.degC)  # reference temperature in degree celsius
    parameters["mesh_density"] = 2 * ureg("")
    parameters["mesh_density_min"] = 2 * ureg("")

    # simulation time
    full_time = 60 * 60 * 1 * ureg("s")  # simulation time
    time_step = 60 * 20 * ureg("s")  # timestep

    # run simulation
    data = demonstrator_beam(full_time, time_step, parameters)
    print(data)
