import warnings

import dolfin as df
import fenics_concrete
import matplotlib.pyplot as plt
import numpy as np
import pint
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning

from lebedigital.unit_registry import ureg

# df.parameters["form_compiler"]["representation"] = "quadrature"
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)


def three_point_bending_example(parameters, pv_output=False):
    """Example of a linear elastic three point bending test

    Parameters
    ----------
        parameters : dict
            - E : float, pint stress unit
                Youngs modulus
            - nu : float, pint unitless
                Poisson's ratio
            Optional
            - height: float, pint length unit
                Height of the beam, default: 300 mm
            - width: float, pint length unit
                Width of the beam, default 150 mm
            - length: float, pint length unit
                Length of the beam, default 2000 mm
            - displacement: float, pint length unit
                Displacement in center of beam, default -10 mm

    Returns
    -------
        stress_in_x : float, pint stress unit in N/mm^2
            Stress in x direction in the center at the bottom, where the maximum stress is expected
    """
    # setting up the simulation parameters
    p = fenics_concrete.Parameters()  # using the current default values

    # default values for the simulation
    p["log_level"] = "WARNING"
    p["dim"] = 3
    p["mesh_density"] = 4  # number of elements in vertical direction

    p["height"] = 300 * ureg("mm")
    p["length"] = 2000 * ureg("mm")
    p["width"] = 150 * ureg("mm")

    # displacement load in the center of the beam
    p["displacement"] = -10 * ureg("mm")

    # will overwrite the default values if given in input
    parameters = p + parameters

    # convert to correct units/check for required inputs
    parameters["height"].ito("mm")
    parameters["length"].ito("mm")
    parameters["width"].ito("mm")
    parameters["displacement"].ito("mm")
    parameters["E"].ito("N/mm^2")
    parameters["nu"].ito("")

    # ... now remove all units, otherwise problems with simulation :(, will change in the future...
    for key in parameters.keys():
        # check if parameter has a unit, then remove it
        if isinstance(parameters[key], type(1 * ureg(""))):
            parameters[key] = parameters[key].magnitude

    # setting up the problem
    experiment = fenics_concrete.ConcreteBeamExperiment(parameters)
    problem = fenics_concrete.LinearElasticity(experiment, parameters)

    problem.experiment.apply_displ_load(parameters["displacement"])

    stress_sensor = fenics_concrete.sensors.StressSensor(df.Point(parameters.length / 2, parameters.width / 2, 0))
    problem.add_sensor(stress_sensor)

    problem.solve(t=0)  # solving this

    if pv_output:
        problem.pv_plot(t=0)

    # results for last (only) load step
    stress_tensor = problem.sensors[stress_sensor.name].data[-1]
    stress_in_x = stress_tensor[0]

    return stress_in_x * ureg("N/mm^2")


if __name__ == "__main__":
    # example of how to use this function
    # defining the material parameters
    parameters = fenics_concrete.Parameters()  # using the current default values

    # input values for the material
    parameters["E"] = 30000 * ureg("N/mm^2")
    parameters["nu"] = 0.2 * ureg("")

    # resulting stress in x direction in the bottom center of the beam
    print("Stress in x:", three_point_bending_example(parameters))
