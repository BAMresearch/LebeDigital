import fenics_concrete
# import probeye
from probeye.definition.forward_model import ForwardModelBase
from probeye.definition.sensor import Sensor


class LinearElasticityCylinder(ForwardModelBase):
    """Probeye forward model for a compression test on a linear elastic cylinder"""

    def interface(self):
        """Definition of the variable parameter, the input and output sensors
        E                 : Young's modulus in N/mm²
        nu                : Poisson's ratio
        height            : height of the cylinder in mm
        radius            : radius of the experimental cylinder in mm
        displacement_list : list with loading conditions in mm
        force_list        : resulting computed forces
        """

        self.parameters = ['E']
        self.input_sensors = [Sensor("nu"),
                              Sensor("height"),
                              Sensor("radius"),
                              Sensor("displacement_list")]
        self.output_sensors = [Sensor('force_list', std_model="sigma")]

    def response(self, inp: dict) -> dict:
        """Setup of the FEM problem

        Parameters
        ----------
            inp : dictionary
                A dictionary with all input values as defined in definition()

        Returns
        -------
            dictionary
                Returns "force_list" as output sensor
        """
        # this method *must* be provided by the user
        parameters = fenics_concrete.Parameters()
        # input parameters
        parameters['E'] = inp["E"]
        parameters['nu'] = inp["nu"]
        parameters['height'] = inp["height"]
        parameters['radius'] = inp["radius"]
        # problem parameters
        parameters['mesh_density'] = 6
        parameters['log_level'] = 'WARNING'
        parameters['bc_setting'] = 'fixed'
        parameters['dim'] = 3
        # setup simulation
        experiment = fenics_concrete.ConcreteCylinderExperiment(parameters)
        problem = fenics_concrete.LinearElasticity(experiment, parameters)
        # setup sensor
        sensor = fenics_concrete.sensors.ReactionForceSensorBottom()
        problem.add_sensor(sensor)

        # loop over all measured displacements
        for displacement in inp["displacement_list"]:
            # apply displacement
            problem.experiment.apply_displ_load(displacement)
            # solve problem
            problem.solve()  # solving this

        # return a list with the measured reaction forces
        force_list = problem.sensors.ReactionForceSensorBottom.data

        return {'force_list': force_list}