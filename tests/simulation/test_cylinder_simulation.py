import os
from pathlib import Path

import fenicsxconcrete
import numpy as np
import pytest
from fenicsxconcrete.experimental_setup.compression_cylinder import CompressionCylinder
from fenicsxconcrete.finite_element_problem.linear_elasticity import LinearElasticity
from fenicsxconcrete.sensor_definition.reaction_force_sensor import ReactionForceSensor
from fenicsxconcrete.util.unit_registry import ureg


def test_cylinder_simulation():
    """Testing the cylinder simulation
    This test is checking if the conda package is correctly installed"""
    parameters = {}  # using the current default values

    parameters["E"] = 3000 * ureg("MPa")
    parameters["nu"] = 0.2 * ureg("")
    parameters["radius"] = 75 * ureg("mm")
    parameters["height"] = 300 * ureg("mm")
    parameters["dim"] = 3 * ureg("")

    parameters["log_level"] = "WARNING" * ureg("")
    parameters["bc_setting"] = "free" * ureg("")
    parameters["mesh_density"] = 6 * ureg("")

    displacement = -3 * ureg("mm")

    experiment = CompressionCylinder(parameters)

    problem = LinearElasticity(experiment, parameters)
    sensor = ReactionForceSensor()

    problem.add_sensor(sensor)

    problem.experiment.apply_displ_load(displacement)

    problem.solve()  # solving this

    # remove temp mesh files that are generated during mesh generation
    directory = Path("mesh")
    files = ["cylinder.h5", "cylinder.msh", "cylinder.xdmf"]

    # delete files
    for file in files:
        path = directory / file
        if path.is_file():
            os.remove(path)

    # delete directory
    os.rmdir(directory)

    # last measurement
    measured_force = problem.sensors[sensor.name].data[-1]

    # exact solution for free bc
    exact_solution = np.pi * parameters["radius"] ** 2 * parameters["E"] * displacement / parameters["height"]

    # scaling factor due to an error resulting from the discretization
    rel_discretization_error = 1.012722

    assert measured_force == pytest.approx(rel_discretization_error * exact_solution)
