import dolfin as df

from usecases.Concrete.simulation_model.material_problems.template_material import MaterialProblem
from usecases.Concrete.simulation_model.helpers import Parameters
from usecases.Concrete.simulation_model import experimental_setups

import warnings
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning

df.parameters["form_compiler"]["representation"] = "quadrature"
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)


# full concrete model, including hydration-temperate and mechanics, including calls to solve etc.
class LinearElasticity(MaterialProblem):
    def __init__(self, experiment=None, parameters=None, pv_name='pv_output_linear_elasticity'):
        # generate "dummy" experiment when none is passed
        if experiment is None:
            experiment = experimental_setups.MinimalCubeExperiment(parameters)

        super().__init__(experiment, parameters, pv_name)

    def setup(self):
        # setup initial temperature material parameters
        default_p = Parameters()
        # polynomial degree
        default_p['degree'] = 2  # default boundary setting

        # parameters for mechanics problem
        default_p['E'] = None  # Young's Modulus
        default_p['nu'] = None  # Poisson's Ratio
        default_p['mu'] = None
        default_p['lmbda'] = None

        self.p = default_p + self.p

        # expecting E and nu to compute mu and lambda, however you can directly supply mu and lambda
        # compute material parameters
        if self.p.mu is None or self.p.lmbda is None:
            assert self.p.E is not None and self.p.nu is not None
            self.p.mu = self.p.E / (2.0 * (1.0 + self.p.nu))
            self.p.lmbda = self.p.E * self.p.nu / ((1.0 + self.p.nu) * (1.0 - 2.0 * self.p.nu))

        # initialize possible paraview output
        self.pv_file = df.XDMFFile(self.pv_name + '.xdmf')
        self.pv_file.parameters["flush_output"] = True
        self.pv_file.parameters["functions_share_mesh"] = True

        # define function space ets.
        self.V = df.VectorFunctionSpace(self.experiment.mesh, "Lagrange", self.p.degree)  # 2 for quadratic elements

        self.residual = None  # initialize residual
        # Define variational problem
        u_trial = df.TrialFunction(self.V)
        v = df.TestFunction(self.V)
        self.a = df.inner(self.sigma(u_trial), df.grad(v)) * df.dx

        if self.p.dim == 2:
            f = df.Constant((0, 0))
        elif self.p.dim == 3:
            f = df.Constant((0, 0, 0))
        else:
            raise Exception(f'wrong dimension {self.p.dim} for problem setup')

        self.L = df.inner(f, v) * df.dx

        # boundary conditions only after function space
        self.bcs = self.experiment.create_displ_bcs(self.V)

        # displacement field
        self.displacement = df.Function(self.V)
        # TODO better names!!!!
        self.visu_space_T = df.TensorFunctionSpace(self.experiment.mesh, "Lagrange", self.p.degree)

    # Stress computation for linear elastic problem
    def sigma(self, v):
        # v is the displacement field
        return 2.0 * self.p.mu * df.sym(df.grad(v)) + self.p.lmbda * df.tr(df.sym(df.grad(v))) * df.Identity(len(v))

    def solve(self, t=1.0):
        # time in this example only relevant for the naming of the paraview steps and the sensor output
        # solve
        df.solve(self.a == self.L, self.displacement, self.bcs)

        # TODO make some switch in sensor definition to trigger this...
        self.compute_residual()

        # get sensor data
        for sensor_name in self.sensors:
            # go through all sensors and measure
            self.sensors[sensor_name].measure(self, t)

    def compute_residual(self):
        # compute reaction forces
        self.residual = df.action(self.a, self.displacement) - self.L

    def pv_plot(self, t=0):
        # paraview output

        # displacement plot
        u_plot = df.project(self.displacement, self.V)
        u_plot.rename("Displacement", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(u_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

        # stress plot
        sigma_plot = df.project(self.sigma(self.displacement), self.visu_space_T)
        sigma_plot.rename("Stress", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(sigma_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
