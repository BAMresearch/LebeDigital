import dolfin as df
import numpy as np


from usecases.Concrete.simulation_model.helpers import Parameters

class Experiment:
    def __init__(self, parameters = None):
        # setup of paramter field
        self.p = Parameters()
        # constants
        self.p['zero_C'] = 273.15  # to convert celcius to kelvin input to

        self.p = self.p + parameters

        self.setup()

    def setup(self):
        raise NotImplementedError()


    # define some common boundary conditions
    def boundary_full(self):
        def bc_full(x, on_boundary):
            return on_boundary
        return bc_full


    def boundary_empty(self):
        def bc_empty(x, on_boundary):
            return None
        return bc_empty


    def boundary_left(self, end = None):
        #print(self.p.dim)
        # Left and right are defined as x (x[0]) in 1D, 2D and 3D
        # minimum vales aka "left" boundary
        left = np.amin(self.mesh.coordinates()[:,0])
        def bc_left(x, on_boundary):
            return on_boundary and df.near(x[0], left)
        return bc_left


    def boundary_right(self):
        # Left and right are defined as x (x[0]) in 1D, 2D and 3D
        # max vales aka "right" boundary
        right = np.amax(self.mesh.coordinates()[:,0])

        def bc_right(x, on_boundary):
            return on_boundary and df.near(x[0], right)
        return bc_right


    def boundary_bottom(self, end = None):
        # the arg "end" is mainly to make a more interesting testing scenario
        # top and bottom are defined as y (x[1]) in 2D and as z (x[2]) in 3D
        if self.p.dim == 2:
            dir_id = 1
        elif self.p.dim == 3:
            dir_id = 2
        else:
            raise Exception('Dimension not defined')

        # minimum vales aka "bottom" boundary
        bottom = np.amin(self.mesh.coordinates()[:,dir_id])
        if end == None:
            end = np.amax(self.mesh.coordinates()[:,0])

        def bc_bottom(x, on_boundary):
            return on_boundary and df.near(x[dir_id], bottom) and x[0] <= end
        return bc_bottom


    def boundary_top(self):
        # top and bottom are defined as y (x[1]) in 2D and as z (x[2]) in 3D
        if self.p.dim == 2:
            dir_id = 1
        elif self.p.dim == 3:
            dir_id = 2
        else:
            raise Exception('Dimension not defined')

            # minimum vales aka "bottom" boundary
        top = np.amax(self.mesh.coordinates()[:, dir_id])

        def bc_top(x, on_boundary):
            return on_boundary and df.near(x[dir_id], top)

        return bc_top


    def boundary_front(self):
        # front and back are not defined in in 2D and as y (x[1]) in 3D
        if self.p.dim == 2:
            bc = self.boundary_empty()

        elif self.p.dim == 3:
            # minimum vales aka "front" boundary
            front = np.amin(self.mesh.coordinates()[:, 1])
            def bc_front(x, on_boundary):
                return on_boundary and df.near(x[1], front)

            bc = bc_front
        else:
            raise Exception('Dimension not defined')

        return bc


    def boundary_back(self):
        # front and back are not defined in in 2D and as y (x[1]) in 3D
        if self.p.dim == 2:
            bc = self.boundary_empty()

        elif self.p.dim == 3:
            # max vales aka "back" boundary
            back = np.amax(self.mesh.coordinates()[:, 1])

            def bc_back(x, on_boundary):
                return on_boundary and df.near(x[1], back)

            bc = bc_back
        else:
            raise Exception('Dimension not defined')

        return bc
