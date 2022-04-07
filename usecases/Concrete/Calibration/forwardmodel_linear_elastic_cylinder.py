import usecases.Concrete.simulation_model as Simulation
#import probeye
from probeye.definition.forward_model import ForwardModelBase
from probeye.definition.sensor import Sensor

class LinearElasticityCylinder(ForwardModelBase):
    def definition(self):
        self.parameters = ['E']
        self.input_sensors = [Sensor("nu"),
                              Sensor("height"),
                              Sensor("radius"),
                              Sensor("displacement_list")]
        self.output_sensors = [Sensor('force_list')]

    def response(self, inp: dict) -> dict:
        # this method *must* be provided by the user
        parameters = Simulation.Parameters()
        # input paramters
        parameters['E'] = inp["E"]
        parameters['nu'] = inp["nu"]
        parameters['height'] = inp["height"]
        parameters['radius'] = inp["radius"]
        # problem paramters
        parameters['mesh_density'] = 6
        parameters['log_level'] = 'WARNING'
        parameters['bc_setting'] = 'fixed'
        parameters['dim'] = 3
        # setup simulation
        experiment = Simulation.Setups.ConcreteCylinderExperiment(parameters)
        problem = Simulation.Models.LinearElasticity(experiment, parameters)
        # setup sensor
        sensor = Simulation.sensors.ReactionForceSensor()
        problem.add_sensor(sensor)

        # loop over all measured displacements
        for displacement in inp["displacement_list"]:
            # apply dispacement
            problem.experiment.apply_displ_load(displacement)
            # solve problem
            problem.solve()  # solving this

        # return a list with the measured recation forces
        force_list = problem.sensors.ReactionForceSensor.data

        return {'force_list': force_list}