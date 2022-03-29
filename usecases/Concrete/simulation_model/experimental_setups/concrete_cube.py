
from usecases.Concrete.simulation_model.experimental_setups.template_experiment import Experiment
from usecases.Concrete.simulation_model.helpers import Parameters
import dolfin as df

class ConcreteCubeExperiment(Experiment):
    def __init__(self, parameters=None):
        # initialize a set of "basic paramters" (for now...)
        p = Parameters()
        # boundary values...
        p['T_0'] = 20  # inital concrete temperature
        p['T_bc1'] = 30  # temperature boundary value 1
        p['T_bc2'] = 50  # temperature boundary value 2
        p['T_bc3'] = 20  # temperature boundary value 3
        p['bc_setting'] = 'full'  # default boundary setting
        p['dim'] = 3  # default boundary setting
        p['mesh_density'] = 10  # default boundary setting
        p['mesh_setting'] = 'left/right'  # default boundary setting
        p = p + parameters
        super().__init__(p)


    def setup(self, bc='full'):
        self.bc = bc  # different boundary settings
        # elements per spacial direction
        n = self.p.mesh_density
        if self.p.dim == 2:
            self.mesh = df.UnitSquareMesh(n, n, self.p.mesh_setting)
        elif self.p.dim == 3:
            self.mesh = df.UnitCubeMesh(n, n, n)
        else:
            print(f'wrong dimension {self.p.dim} for problem setup')
            exit()

        # self.p['T_0'] = 20 # inital concrete temperature
        # self.p['T_bc1'] = 10 # temperature boundary value 1
        # self.p['T_bc2'] = 50 # temperature boundary value 2
        # self.p['T_bc3'] = 10 # temperature boundary value 3

    def create_temp_bcs(self, V):
        # Temperature boundary conditions
        T_bc1 = df.Expression('t_boundary', t_boundary=self.p.T_bc1 + self.p.zero_C, degree=0)
        T_bc2 = df.Expression('t_boundary', t_boundary=self.p.T_bc2 + self.p.zero_C, degree=0)
        T_bc3 = df.Expression('t_boundary', t_boundary=self.p.T_bc3 + self.p.zero_C, degree=0)

        temp_bcs = []

        if self.p.bc_setting == 'full':
            # bc.append(DirichletBC(temperature_problem.V, T_bc, full_boundary))
            temp_bcs.append(df.DirichletBC(V, T_bc1, self.boundary_full()))
        elif self.p.bc_setting == 'test-setup':
            # bc.append(DirichletBC(temperature_problem.V, T_bc, full_boundary))
            temp_bcs.append(df.DirichletBC(V, T_bc1, self.boundary_left()))
            temp_bcs.append(df.DirichletBC(V, T_bc1, self.boundary_bottom(0.5)))
            temp_bcs.append(df.DirichletBC(V, T_bc2, self.boundary_right()))
        else:
            raise Exception(
                f'parameter[\'bc_setting\'] = {self.bc_setting} is not implemented as temperature boundary.')

        return temp_bcs

    def create_displ_bcs(self, V):
        # define displacement boundary
        displ_bcs = []

        if self.p.dim == 2:
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0)), self.boundary_bottom()))
        elif self.p.dim == 3:
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0, 0)),  self.boundary_bottom()))

        return displ_bcs