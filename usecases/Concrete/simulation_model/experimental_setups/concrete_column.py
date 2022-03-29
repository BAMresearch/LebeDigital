
from usecases.Concrete.simulation_model.experimental_setups.template_experiment import Experiment
from usecases.Concrete.simulation_model.helpers import Parameters
import numpy as np
import dolfin as df

class ConcreteColumnExperiment(Experiment):
    def __init__(self, parameters=None):
        # initialize a set of "basic paramters" (for now...)
        p = Parameters()
        # boundary values...
        p['T_0'] = 20  # inital concrete temperature
        p['T_bc1'] = 40  # temperature boundary value 1
        p['bc_setting'] = 'full'  # default boundary setting
        p['dim'] = 3  # default boundary setting
        p['mesh_density'] = 4  # default boundary setting
        p['mesh_density_min'] = 1
        p['mesh_setting'] = 'left/right'  # default boundary setting

        p['width'] = 0.5  # length of pillar in m
        p['height'] = 4 # width (square crossection)

        p = p + parameters
        super().__init__(p)






    def setup(self):
        # minimale number of elements
        md_width = np.amax([self.p.mesh_density_min,int(self.p.mesh_density * self.p.width)])
        md_height = np.amax([self.p.mesh_density_min,int(self.p.mesh_density * self.p.height)])

        if self.p.dim == 2:
             self.mesh = df.RectangleMesh(df.Point(0., 0.), df.Point(self.p.width, self.p.height),
                                          md_width, md_height, diagonal='right')

        elif self.p.dim == 3:
            self.mesh = df.BoxMesh(df.Point(0, 0, 0), df.Point(self.p.width, self.p.width, self.p.height),
                                   md_width, md_width, md_height)
        else:
            raise Exception(f'wrong dimension {self.p.dim} for problem setup')


    def create_temp_bcs(self, V):
        # Temperature boundary conditions
        T_bc1 = df.Expression('t_boundary', t_boundary=self.p.T_bc1 + self.p.zero_C, degree=0)

        temp_bcs = []

        if self.p.bc_setting == 'full':
             temp_bcs.append(df.DirichletBC(V, T_bc1, self.boundary_full()))
        else:
             raise Exception(
                 f'parameter[\'bc_setting\'] = {self.p.bc_settings} is not implemented as temperature boundary.')

        return temp_bcs


    def create_displ_bcs(self, V):
        # define displacement boundary
        displ_bcs = []

        if self.p.dim == 2:
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0)), self.boundary_bottom()))
        elif self.p.dim == 3:
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0, 0)),  self.boundary_bottom()))

        return displ_bcs