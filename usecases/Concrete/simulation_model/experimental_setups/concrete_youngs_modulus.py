
from usecases.Concrete.simulation_model.experimental_setups.template_experiment import Experiment
from usecases.Concrete.simulation_model.helpers import Parameters
import numpy as np
import dolfin as df
import mshr

class ConcreteYoungsModulusExperiment(Experiment):
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

        p['radius'] = 75  # length of pillar in m
        p['height'] = 100 # width (square crossection)

        p = p + parameters
        super().__init__(p)






    def setup(self):

        # if self.p.dim == 2:
        #      self.mesh = df.RectangleMesh(df.Point(0., 0.), df.Point(self.p.width, self.p.height),
        #                                   md_width, md_height, diagonal='right')

        if self.p.dim == 3:

            # The mesh geometry
            # Cylinder ( center bottom, center top, radius bottom, radius top )
            cylinder_geometry = mshr.Cylinder(df.Point(0, 0, 0), df.Point(0, 0, self.p.height), self.p.radius, self.p.radius)

            # mesh ( geometry , mesh density )
            self.mesh = mshr.generate_mesh(cylinder_geometry, self.p.mesh_density)



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
        # todo define this somewhere else
        # self.top_displacement.assign(Constant(new value))
        self.top_displacement = df.Constant(50.0)

        if self.p.dim == 2:
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0)), self.boundary_bottom()))
        elif self.p.dim == 3:
            displ_bcs.append(df.DirichletBC(V.sub(2), self.top_displacement, self.boundary_top()))  # apply displacement
            displ_bcs.append(df.DirichletBC(V.sub(0), 0, self.boundary_top()))
            displ_bcs.append(df.DirichletBC(V.sub(1), 0, self.boundary_top()))
            #displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0, )),  self.boundary_top()))
            displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0, 0)),  self.boundary_bottom()))

        return displ_bcs