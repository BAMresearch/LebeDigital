

from usecases.Concrete.simulation_model.experimental_setups.template_experiment import Experiment
from usecases.Concrete.simulation_model.helpers import Parameters
import dolfin as df


class ConcreteBeamExperiment(Experiment):
    def __init__(self, parameters = None):
        p = Parameters()
        # boundary values...
        p['T_0'] = 20  # inital concrete temperature
        p['T_bc1'] = 30  # temperature boundary value 1
        p['T_bc2'] = 50  # temperature boundary value 2
        p['T_bc3'] = 20  # temperature boundary value 3
        p['length'] = 5 # m length
        p['height'] = 1 # height
        p['width'] = 0.8 # width
        p['mesh_density'] = 4  # default boundary setting
        p['bc_setting'] = 'full' # default boundary setting
        p = p + parameters
        super().__init__(p)

    def setup(self):
        # elements per spacial direction

        if self.p.dim == 2:
             self.mesh = df.RectangleMesh(df.Point(0., 0.), df.Point(self.p.length, self.p.height),
                                          int(self.p.mesh_density * self.p.length),
                                          int(self.p.mesh_density * self.p.height), diagonal='right')
        elif self.p.dim == 3:
            self.mesh = df.BoxMesh(df.Point(0, 0, 0), df.Point(self.p.length, self.p.width, self.p.height),
                                   int(self.p.length * self.p.mesh_density),
                                   int(self.p.width * self.p.mesh_density),
                                   int(self.p.height * self.p.mesh_density))
        else:
            raise Exception(f'wrong dimension {self.p.dim} for problem setup')


    def create_temp_bcs(self,V):

        # Temperature boundary conditions
        T_bc1 = df.Expression('t_boundary', t_boundary=self.p.T_bc1+self.p.zero_C, degree=0)
        T_bc2 = df.Expression('t_boundary', t_boundary=self.p.T_bc2+self.p.zero_C, degree=0)
        T_bc3 = df.Expression('t_boundary', t_boundary=self.p.T_bc3+self.p.zero_C, degree=0)

        temp_bcs = []

        if self.p.bc_setting == 'full':
            # bc.append(DirichletBC(temperature_problem.V, T_bc, full_boundary))
            temp_bcs.append(df.DirichletBC(V, T_bc1, self.boundary_full()))
        elif self.p.bc_setting == 'left-right':
            # bc.append(DirichletBC(temperature_problem.V, T_bc, full_boundary))
            temp_bcs.append(df.DirichletBC(V, T_bc2, self.boundary_left()))
            temp_bcs.append(df.DirichletBC(V, T_bc3, self.boundary_right()))
        else:
            raise Exception(f'parameter[\'bc_setting\'] = {self.p.bc_setting} is not implemented as temperature boundary.')

        return temp_bcs


    def create_displ_bcs(self,V):
        if self.p.dim == 2:
            dir_id = 1
            fixed_bc = df.Constant((0, 0))
        elif self.p.dim == 3:
            dir_id = 2
            fixed_bc = df.Constant((0, 0, 0))

        # define surfaces, full, left, right, bottom, top, none
        def left_support(x, on_boundary):
            return df.near(x[0], 0) and df.near(x[dir_id], 0)
        def right_support(x, on_boundary):
            return df.near(x[0], self.p.l) and df.near(x[dir_id], 0)

        # define displacement boundary
        displ_bcs = []
        displ_bcs.append(df.DirichletBC(V, fixed_bc, left_support, method='pointwise'))
        displ_bcs.append(df.DirichletBC(V.sub(dir_id), df.Constant(0), right_support, method='pointwise'))

        return displ_bcs
        
  
