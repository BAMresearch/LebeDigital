import dolfin as df
import numpy as np
import scipy.optimize


from usecases.Concrete.simulation_model.material_problems.template_material import MaterialProblem

from usecases.Concrete.simulation_model.helpers import Parameters
from usecases.Concrete.simulation_model.helpers import set_q
from usecases.Concrete.simulation_model.helpers import LocalProjector
from usecases.Concrete.simulation_model import experimental_setups


import warnings
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning

df.parameters["form_compiler"]["representation"] = "quadrature"
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)


# full concrete model, including hydration-temperate and mechanics, including calls to solve etc.
class LinearElasticity(MaterialProblem):
    def __init__(self, experiment=None, parameters=None, pv_name='pv_output_linear_elasticity'):
        # generate "dummy" experiement when none is passed
        if experiment == None:
            experiment = experimental_setups.get_experiment('MinimalCube', parameters)

        super().__init__(experiment, parameters, pv_name)

    def setup(self):
        # setup initial temperatre material paramters
        default_p = Parameters()
        # polinomial degree
        default_p['degree'] = 2  # default boundary setting

        ### paramters for mechanics problem
        default_p['E'] = None # Youngs Modulus
        default_p['nu'] = None  # Poissons Ratio
        default_p['mu'] = None
        default_p['lmbda'] = None

        self.p = default_p + self.p

        # expecting E and nu to compute mu and lambda, however you can directly supply mu and lambda
        # compute material paramters
        if self.p.mu == None or self.p.lmbda == None:
            assert self.p.E != None and self.p.nu != None
            self.p.mu =  self.p.E / (2.0 * (1.0 + self.p.nu))
            self.p.lmbda = self.p.E * self.p.nu / ((1.0 + self.p.nu) * (1.0 - 2.0 * self.p.nu))

        # initialize possible paraview output
        self.pv_file = df.XDMFFile(self.pv_name + '.xdmf')
        self.pv_file.parameters["flush_output"] = True
        self.pv_file.parameters["functions_share_mesh"] = True

        # define function space ets.
        self.V = df.VectorFunctionSpace(self.experiment.mesh, "Lagrange", self.p.degree)  # 2 for quadratic elements

        # Define variational problem
        u_trial = df.TrialFunction(self.V)
        v = df.TestFunction(self.V)
        self.a = df.inner(self.sigma(u_trial), df.grad(v)) * df.dx
        f = df.Constant((0, 0, 0))
        self.L = df.inner(f, v) * df.dx

        # bounary conditions only after function space
        self.bcs = self.experiment.create_displ_bcs(self.V)

        # displacement field
        self.u = df.Function(self.V)
        # TODO better namess!!!!
        self.visu_space_T = df.TensorFunctionSpace(self.experiment.mesh, "Lagrange", self.p.degree)

    # Stress computation for linear elastic problem
    def sigma(self, v):
        # v is the displacement field
        return 2.0 * self.p.mu * df.sym(df.grad(v)) + self.p.lmbda * df.tr(df.sym(df.grad(v))) * df.Identity(len(v))

    def solve(self, t=1.0):
        # TODO: apply the displecement load...
        #       - how do I define "steps"?,

        # solve
        df.solve(self.a == self.L, self.u, self.bcs)

        # TODO: do things for residual computation...
        #       - implement force sensor

        # get sensor data
        for sensor_name in self.sensors:
            # go through all sensors and measure
            self.sensors[sensor_name].measure(self, t)

    def pv_plot(self, t=0):
        # paraview output

        # displacement plot
        u_plot = df.project(self.u, self.V)
        u_plot.rename("Displacement", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(u_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

        # stress plot
        sigma_plot = df.project(self.sigma(self.u), self.visu_space_T)
        sigma_plot.rename("Stress", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(sigma_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

    def set_timestep(self, dt):
        # TODO maybe implement a set bc_displacement???
        #      - or is this a function of the "experiment"???
        #self.temperature_problem.set_timestep(dt)
        pass




