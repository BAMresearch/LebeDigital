from dolfin import *
from mshr import *
import numpy as np

#assuming stress/strain data

class FenicsYoungsModulusTestSimulation():
    def __init__(self):
        pass

    def __call__(self, inp):
        """
                Simple simulation representing a Young's modulus test for concrete.
                Based on an applied strain, the stress in loading direction in the center is computed.

                Input
                ----------
                inp : Input library
                        inp['E'] Young's modulus for simulation
                        inp['nu'] Poisson's ratio for simulation
                        inp['height'] Height of test cylinder in mm
                        optional:
                            inp['area'] crossection of test cylinder in mm^2
                            inp['radius'] Radius of test cylinder in mm
                            inp['diameter'] Diameter of test cylinder in mm
                        optional (if both are passed displacements are used):
                          inp['displ'] A list with all displ to be applied
                          inp['strains'] A list with all strains to be applied

                Output
                ----------
                response: Output library
                            response['force'] A list with computed reaction forces
                            response['stress'] A list with computed stresses from reaction forces

        """
        # TODO:
        # - how to handle units?
        # - include a check and errors for missing parameters?
        # - compute real discretized surface area

        # set boundary option
        # true : no horzontal displacement at boundary - more realistic option
        # false : horizontal displacement possile - to compare to analytical result (homogeneous stress)
        constrained_bc = True

        # Define and check inputs
        E = inp['E']                    # Youngs modulus
        nu = inp['nu']                  # Poisson's ratio, this can't be calibrated based on this test
        mesh_density = 10                # Parameter for mesh generation
        height = inp['height']          # Height of test cylinder in mm

        # compute radius, this allows area as input to check the resulting force/stress more intuitively
        if 'area' in inp:
            radius = (inp['area']/np.pi)**0.5
        elif 'radius' in inp:
            radius = inp['radius']
        elif 'diameter' in inp:
            radius = inp['diameter']/2
        else:
            # fancy error message
            pass

        # check loading input
        if 'displ' in inp:
            input_displ = inp['displ']  # List of displacements to be applied
        elif 'strain' in inp:
            # compute displacement from strain list
            input_displ= np.array(inp['strain']) * height
        else:
            # fancy error message
            pass

        # Elasticity parameters
        mu = E / (2.0 * (1.0 + nu))
        lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

        # The mesh geometry
        # Cylinder ( center bottom, center top, radius bottom, radius top )
        cylinder_geometry = Cylinder(Point(0, 0, 0), Point(0, 0, height), radius, radius)

        # mesh ( geometry , mesh density )
        mesh = generate_mesh(cylinder_geometry, mesh_density)

        # Create function space
        V = VectorFunctionSpace(mesh, "Lagrange", 2)  # 2 for quadratic elements

        # define surfaces
        def bottom_surface(x, on_boundary):
            return on_boundary and near(x[2], 0)

        def top_surface(x, on_boundary):
            return on_boundary and near(x[2], height)

        def origin_node(x, on_boundary):
            return on_boundary and near(x[0], 0) and near(x[1], 0) and near(x[2], 0)

        def side_node(x, on_boundary):
            return on_boundary and near(x[1], radius) and near(x[2], 0)

        # loop over all strains
        computed_stress = []
        computed_force = []
        for top_displacement in input_displ:
            # negative displacement values are compression

            # Set up boundary conditions
            if constrained_bc:
                # dirichlet boundary at top and bottom surface, displacement applied at top
                bc1 = DirichletBC(V.sub(2), top_displacement, top_surface)  # apply displacement
                bc2 = DirichletBC(V.sub(0), 0, top_surface)
                bc3 = DirichletBC(V.sub(1), 0, top_surface)
                bc4 = DirichletBC(V, Constant((0, 0, 0)), bottom_surface)

                # combine all applied boundary conditions
                bcs = [bc1, bc2, bc3, bc4]
            else:
                # dirichlet boundary at top and bottom surface in loading direction only, displacement applied at top
                bc1 = DirichletBC(V.sub(2), top_displacement, top_surface)  # apply displacement
                bc2 = DirichletBC(V.sub(2), Constant(0), bottom_surface)
                bc3 = DirichletBC(V,Constant((0,0,0)),origin_node,method="pointwise") # no cylinder displacement
                bc4 = DirichletBC(V.sub(0),Constant(0),side_node,method="pointwise") # no roation

                # combine all applied boundary conditions
                bcs = [bc1, bc2, bc3, bc4]

            # Stress computation for linear elastic problem
            def sigma(v):
                return 2.0 * mu * sym(grad(v)) + lmbda * tr(sym(grad(v))) * Identity(len(v))

            # Define variational problem
            u = TrialFunction(V)
            v = TestFunction(V)
            a = inner(sigma(u), grad(v)) * dx
            f = Constant((0, 0, 0))
            L = inner(f, v) * dx

            # solve
            u = Function(V)
            solve(a == L, u, bcs)

            # compute reaction forces
            residual = action(a,u) - L
            v_reac = Function(V)
            bc_z = DirichletBC(V.sub(2), Constant(1.), bottom_surface)
            bc_z.apply(v_reac.vector())

            # get forces
            computed_force.append(-assemble(action(residual, v_reac)))

            # compute stress by dividing reaction forces by (approximeted) area
            computed_stress.append(computed_force[-1]/(np.pi*radius**2))

        response = {'stress' : computed_stress,
                     'force' : computed_force}

        return response


#generate simulation object
simulation = FenicsYoungsModulusTestSimulation()

# set input
input = {'E' : 3000,
         'nu' : 0.2,
         'height' : 100,
         'diameter' : 150,
         'displ':[0,-1,-2.5]
         }

# run the model
model_answer = simulation(input)
print(model_answer['force'])
print(model_answer['stress'])


