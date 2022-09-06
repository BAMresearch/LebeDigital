import dolfin as df
import fenics_concrete


def three_point_bending_example(E, nu):
    """Example of a linear elastic three point bending test

    Parameters
    ----------
        E : float
            Young's modulus in N/mm²
        nu : float
            Poisson's ratio

    Returns
    -------
        stress_in_x : float
            Stress in x direction in the center at the bottom, where the maximum stress is expected
    """
    # setting up the simulation parameters
    parameters = fenics_concrete.Parameters()  # using the current default values
    # input values for the material
    parameters['E'] = E
    parameters['nu'] = nu
    # definition of the beam and mesh
    parameters['dim'] = 3
    parameters['mesh_density'] = 4  # number of elements in vertical direction
    parameters['height'] = 300  # in mm
    parameters['length'] = 2000  # in mm
    parameters['width'] = 150  # in mm
    parameters['log_level'] = 'WARNING'
    # displacement load in the center of the beam
    displacement = -10  # displacement load in the center of the beam in mm

    # setting up the problem
    experiment = fenics_concrete.ConcreteBeamExperiment(parameters)
    problem = fenics_concrete.LinearElasticity(experiment, parameters)

    # applying the load
    problem.experiment.apply_displ_load(displacement)

    # applying the stress sensor
    stress_sensor = fenics_concrete.sensors.StressSensor(df.Point(parameters.length/2, 0, 0))
    problem.add_sensor(stress_sensor)

    # solving the problem
    problem.solve(t=0)  # solving this

    # results for last (only) load step
    stress_tensor = problem.sensors[stress_sensor.name].data[-1]
    stress_in_x = stress_tensor[0]

    return stress_in_x


# example of how to use this function
# defining the material parameters
# E = 30000  # N/mm²
# nu = 0.2
# stress = three_point_bending_example(E, nu)
# # resulting stress in x direction in the bottom center of the beam
# print(stress)
