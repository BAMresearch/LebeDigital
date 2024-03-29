import fenics_concrete
import pandas as pd
import pytest

from lebedigital.simulation.precast_column import column_simulation
from lebedigital.unit_registry import ureg


def test_column_simulation():
    # setup parameters:
    parameters = fenics_concrete.Parameters()

    # model parameters

    Q_ = ureg.Quantity
    # concrete values
    parameters["density"] = 2350 * ureg("kg/m^3")  # in kg/m^3 density of concrete
    parameters["themal_cond"] = 2.0 * ureg("W/m/K")  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
    parameters["vol_heat_cap"] = 2.4e6 * ureg("J/m^3/K")  # volumetric heat cap J/m3

    # parameters for mechanics problem
    parameters["E"] = 25e9 * ureg("N/m^2")  # Youngs Modulus in Pa
    parameters["nu"] = 0.2 * ureg("")  # Poissons Ratio

    # required parameters for alpha to E mapping
    parameters["alpha_t"] = 0.2 * ureg("")
    parameters["alpha_0"] = 0.05 * ureg("")
    parameters["a_E"] = 0.6 * ureg("")

    # required parameters for alpha to tensile and compressive stiffness mapping
    parameters["fc"] = 30e6 * ureg("N/m^2")  # in Pa
    parameters["a_fc"] = 1.5 * ureg("")
    parameters["ft"] = 4e6 * ureg("N/m^2")  # in Pa
    parameters["a_ft"] = 1.2 * ureg("")

    # temperature setings:
    parameters["T_0"] = Q_(40, ureg.degC)  # initial temperature of concrete
    parameters["T_bc1"] = Q_(20, ureg.degC)  # constant boundary temperature

    # column geometry
    parameters["width"] = 0.5 * ureg("m")  # length of pillar in m
    parameters["height"] = 4 * ureg("m")  # width (square cross-section)

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

    #
    # simulation time
    full_time = 60 * 60 * 1 * ureg("s")  # simulation time
    time_step = 60 * 20 * ureg("s")  # timestep

    # run simulation
    data = column_simulation(full_time, time_step, parameters)

    assert data.time.values.quantity.magnitude == pytest.approx([1200, 2400, 3600])
    assert data.temperature.values.quantity.magnitude == pytest.approx([41.487825, 43.581025, 48.334999])
    assert data["yield"].values.quantity.magnitude == pytest.approx(
        [21.33535680276793, 5.274222603247022, 1.9011303929033334]
    )
