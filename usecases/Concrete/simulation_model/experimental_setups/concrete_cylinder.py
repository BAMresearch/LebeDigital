from usecases.Concrete.simulation_model.experimental_setups.template_experiment import Experiment
from usecases.Concrete.simulation_model.helpers import Parameters
import dolfin as df
import mshr


class ConcreteCylinderExperiment(Experiment):
    def __init__(self, parameters=None):
        # initialize a set of "basic parameters" (for now...)
        p = Parameters()
        # boundary values...
        p['T_0'] = 20  # initial concrete temperature
        p['T_bc1'] = 40  # temperature boundary value 1
        p['bc_setting'] = 'full'  # default boundary setting
        p['dim'] = 3  # default boundary setting
        p['mesh_density'] = 4  # default boundary setting
        p['mesh_density_min'] = 1
        p['mesh_setting'] = 'left/right'  # default boundary setting
        p['bc_setting'] = 'free'

        p['radius'] = 75   # length of pillar in m
        p['height'] = 100  # width (square cross-section)

        p = p + parameters
        super().__init__(p)

        # initialize variable top_displacement
        self.top_displacement = df.Constant(0.0)  # applied via fkt: apply_displ_load(...)

    def setup(self):
        if self.p.dim == 2:
            md_width = int(self.p.radius*2 * self.p.mesh_density)
            md_height = int(self.p.height * self.p.mesh_density)
            self.mesh = df.RectangleMesh(df.Point(0., 0.), df.Point(self.p.radius*2, self.p.height),
                                         md_width, md_height, diagonal='right')
        elif self.p.dim == 3:
            # The mesh geometry
            # Cylinder ( center bottom, center top, radius bottom, radius top )

            # testing mesh generation....

            cylinder_geometry = mshr.Cylinder(df.Point(0, 0, 0), df.Point(0, 0, self.p.height),
                                              self.p.radius, self.p.radius)

            # mesh ( geometry , mesh density )
            self.mesh = mshr.generate_mesh(cylinder_geometry, self.p.mesh_density)

        else:
            raise Exception(f'wrong dimension {self.p.dim} for problem setup')

    def create_displ_bcs(self, V):
        # define displacement boundary
        displ_bcs = []

        def boundary_node(point):
            def bc_node(x, on_boundary):
                return df.near(x[0], point[0]) and df.near(x[1], point[1]) and df.near(x[2], point[2])
            return bc_node

        if self.p.bc_setting == 'fixed':
            if self.p.dim == 2:
                displ_bcs.append(df.DirichletBC(V.sub(1), self.top_displacement, self.boundary_top()))  # displacement
                displ_bcs.append(df.DirichletBC(V.sub(0), 0, self.boundary_top()))
                displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0)), self.boundary_bottom()))

            elif self.p.dim == 3:
                displ_bcs.append(df.DirichletBC(V.sub(2), self.top_displacement, self.boundary_top()))  # displacement
                displ_bcs.append(df.DirichletBC(V.sub(0), 0, self.boundary_top()))
                displ_bcs.append(df.DirichletBC(V.sub(1), 0, self.boundary_top()))
                displ_bcs.append(df.DirichletBC(V, df.Constant((0, 0, 0)),  self.boundary_bottom()))

        elif self.p.bc_setting == 'free':
            if self.p.dim == 3:
                # getting nodes at the bottom of the mesh to apply correct boundary condition to arbitrary cylinder mesh
                mesh_points = self.mesh.coordinates()  # list of all nodal coordinates
                mesh_points = mesh_points[mesh_points[:, 2].argsort()]  # sorting by z coordinate
                i = 0
                while mesh_points[i][2] == 0.0:
                    i += 1
                bottom_points = mesh_points[:i]
                x_min_boundary_point = bottom_points[bottom_points[:, 0].argsort(kind='mergesort')][
                    0]  # sorting by x coordinate
                x_max_boundary_point = bottom_points[bottom_points[:, 0].argsort(kind='mergesort')][
                    -1]  # sorting by x coordinate
                y_boundary_point = bottom_points[bottom_points[:, 1].argsort(kind='mergesort')][
                    0]  # sorting by y coordinate

                displ_bcs.append(df.DirichletBC(V.sub(2), self.top_displacement, self.boundary_top()))  # displacement
                displ_bcs.append(df.DirichletBC(V.sub(2), 0.9, self.boundary_bottom()))
                displ_bcs.append(df.DirichletBC(V.sub(1), 0.0, boundary_node(x_min_boundary_point), method="pointwise"))
                displ_bcs.append(df.DirichletBC(V.sub(1), 0.0, boundary_node(x_max_boundary_point), method="pointwise"))
                displ_bcs.append(df.DirichletBC(V.sub(0), 0.0, boundary_node(y_boundary_point), method="pointwise"))

        else:
            raise Exception(f'Wrong boundary setting: {self.p.bc_setting}, for cylinder setup')

        return displ_bcs

    def apply_displ_load(self, top_displacement):
        self.top_displacement.assign(df.Constant(top_displacement))
