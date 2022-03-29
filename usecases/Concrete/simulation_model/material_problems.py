import dolfin as df
import numpy as np
import scipy.optimize

from concrete_model.helpers import Parameters
from concrete_model.sensors import Sensors
from concrete_model.helpers import set_q
from concrete_model.helpers import LocalProjector
from concrete_model import experimental_setups

from loguru import logger
import logging
import sys

import warnings
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning

df.parameters["form_compiler"]["representation"] = "quadrature"
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)


class MaterialProblem():
    def __init__(self, experiment, parameters=None, pv_name='pv_output_full'):
        self.experiment = experiment
        # setting up paramters
        self.p = Parameters()
        # constants
        # TODO: where to put these?, what about units???
        self.p['zero_C'] = 273.15
        self.p['igc'] = 8.3145  # ideal gas constant [JK −1 mol −1 ]
        self.p['g'] = 9.81  # graviational acceleration in m/s²

        # other "globel" paramters...
        self.p['log_level'] = 'INFO'

        self.p = self.p + self.experiment.p + parameters

        # set log level...
        if self.p.log_level == 'NOTSET':
            df.set_log_level(0)
            logging.getLogger("FFC").setLevel(logging.NOTSET)
            logging.getLogger("UFL").setLevel(logging.NOTSET)
            logger.add(sys.stderr, level="NOTSET")
        elif self.p.log_level == 'DEBUG':
            df.set_log_level(10)
            logging.getLogger("FFC").setLevel(logging.DEBUG)
            logging.getLogger("UFL").setLevel(logging.DEBUG)
            logger.add(sys.stderr, level="DEBUG")
        elif self.p.log_level == 'INFO':
            df.set_log_level(20)
            logging.getLogger("FFC").setLevel(logging.INFO)
            logging.getLogger("UFL").setLevel(logging.INFO)
            logger.add(sys.stderr, level="INFO")
        elif self.p.log_level == 'WARNING':
            df.set_log_level(30)
            logging.getLogger("FFC").setLevel(logging.WARNING)
            logging.getLogger("UFL").setLevel(logging.WARNING)
            logger.add(sys.stderr, level="WARNING")
        elif self.p.log_level == 'ERROR':
            df.set_log_level(40)
            logging.getLogger("FFC").setLevel(logging.ERROR)
            logging.getLogger("UFL").setLevel(logging.ERROR)
            logger.add(sys.stderr, level="ERROR")
        elif self.p.log_level == 'CRITICAL':
            df.set_log_level(50)
            logging.getLogger("FFC").setLevel(logging.CRITICAL)
            logging.getLogger("UFL").setLevel(logging.CRITICAL)
            logger.add(sys.stderr, level="CRITICAL")
        else:
            level = self.p['log_level']
            raise Exception(f'unknown log level {level}')


        self.sensors =  Sensors()  # list to hold attached sensors


        self.pv_name = pv_name

        #setup fields for sensor output, can be defined in model
        self.displacement = None
        self.temperature = None
        self.degree_of_hydration = None
        self.q_degree_of_hydration = None

        # setup the material object to access the function
        self.setup()

    def setup(self):
        # initialization of this specific problem
        raise NotImplementedError()

    def solve(self):
        # define what to do, to solve this problem
        raise NotImplementedError()

    def add_sensor(self, sensor):

        self.sensors[sensor.name] = sensor





# full concrete model, including hydration-temperate and mechanics, including calls to solve etc.
class ConcreteThermoMechanical(MaterialProblem):
    def __init__(self, experiment=None, parameters=None, pv_name='pv_output_concrete-thermo-mechanical'):
        # generate "dummy" experiement when none is passed
        if experiment == None:
            experiment = experimental_setups.get_experiment('MinimalCube', parameters)

        super().__init__(experiment, parameters, pv_name)

    #     # TODO: define global fields here
    #     #       - alpha, V
    #     #       - etc...

    def setup(self):
        # setup initial temperatre material paramters
        default_p = Parameters()
        # Material parameter for concrete model with temperature and hydration
        default_p['density'] = 2350  # in kg/m^3 density of concrete
        default_p['density_binder'] = 1440  # in kg/m^3 density of the binder
        default_p['themal_cond'] = 2.0  # effective thermal conductivity, approx in Wm^-3K^-1, concrete!
        # self.specific_heat_capacity = 9000  # effective specific heat capacity in J kg⁻1 K⁻1
        default_p['vol_heat_cap'] = 2.4e6  # volumetric heat cap J/(m3 K)
        default_p['b_ratio'] = 0.2  # volume percentage of binder
        default_p['Q_pot'] = 500e3  # potential heat per weight of binder in J/kg
        # p['Q_inf'] = self.Q_pot * self.density_binder * self.b_ratio  # potential heat per concrete volume in J/m3
        default_p['B1'] = 2.916E-4  # in 1/s
        default_p['B2'] = 0.0024229  # -
        default_p['eta'] = 5.554  # something about diffusion
        default_p['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c
        default_p['E_act'] = 5653 * self.p.igc  # activation energy in Jmol^-1
        default_p['T_ref'] = 25  # reference temperature in degree celsius
        # setting for temperature adjustment
        # option: 'exponential' and 'off'
        default_p['temp_adjust_law'] = 'exponential'
        # polinomial degree
        default_p['degree'] = 2  # default boundary setting

        ### paramters for mechanics problem
        default_p['E_28'] = 15000000  # Youngs Modulus N/m2 or something... TODO: check units!
        default_p['nu'] = 0.2  # Poissons Ratio

        # required paramters for alpha to E mapping
        default_p['alpha_t'] = 0.2
        default_p['alpha_0'] = 0.05
        default_p['a_E'] = 0.6

        # required paramters for alpha to tensile and compressive stiffness mapping
        default_p['fc_inf'] = 6210000
        default_p['a_fc'] = 1.2
        default_p['ft_inf'] = 467000
        default_p['a_ft'] = 1.0

        self.p = default_p + self.p

        # setting up the two nonlinear problems
        self.temperature_problem = ConcreteTempHydrationModel(self.experiment.mesh, self.p, pv_name=self.pv_name)

        # here I "pass on the parameters from temperature to mechanics problem.."
        self.mechanics_problem = ConcreteMechanicsModel(self.experiment.mesh, self.p, pv_name=self.pv_name)
        # coupling of the output files
        self.mechanics_problem.pv_file = self.temperature_problem.pv_file

        # initialize concrete temperature as given in experimental setup
        self.set_inital_T(self.p.T_0)

        # setting bcs
        self.mechanics_problem.set_bcs(self.experiment.create_displ_bcs(self.mechanics_problem.V))
        self.temperature_problem.set_bcs(self.experiment.create_temp_bcs(self.temperature_problem.V))

        # setting up the solvers
        self.temperature_solver = df.NewtonSolver()
        self.temperature_solver.parameters['absolute_tolerance'] = 1e-9
        self.temperature_solver.parameters['relative_tolerance'] = 1e-8

        self.mechanics_solver = df.NewtonSolver()
        self.mechanics_solver.parameters['absolute_tolerance'] = 1e-9
        self.mechanics_solver.parameters['relative_tolerance'] = 1e-8





    def solve(self, t=1.0):

        # print('Solving: T') # TODO ouput only a certain log level INFO
        self.temperature_solver.solve(self.temperature_problem, self.temperature_problem.T.vector())

        # set current DOH for computation of Young's modulus
        self.mechanics_problem.q_alpha = self.temperature_problem.q_alpha
        # print('Solving: u') # TODO ouput only a certain log level INFO

        # mechanics paroblem is not required for temperature, could crash in frist time steps but then be useful
        try:
            self.mechanics_solver.solve(self.mechanics_problem, self.mechanics_problem.u.vector())
        except Exception as e:
            print(f'Mechanics crashed at time: {t}')
            print(e)

        # history update
        self.temperature_problem.update_history()

        # save fields to global problem for sensor output
        self.displacement = self.mechanics_problem.u
        self.temperature = self.temperature_problem.T
        self.degree_of_hydration = df.project(self.temperature_problem.q_alpha, self.temperature_problem.visu_space, form_compiler_parameters={'quadrature_degree': self.p.degree})
        self.q_degree_of_hydration = self.temperature_problem.q_alpha
        self.q_yield = self.mechanics_problem.q_yield

        # get sensor data
        for sensor_name in self.sensors:
            # go through all sensors and measure
            self.sensors[sensor_name].measure(self, t)

    def pv_plot(self, t=0):
        # calls paraview output for both problems
        self.temperature_problem.pv_plot(t=t)
        self.mechanics_problem.pv_plot(t=t)

    def set_inital_T(self, T):
        self.temperature_problem.set_initial_T(T)

    def set_timestep(self, dt):
        self.temperature_problem.set_timestep(dt)

    def get_heat_of_hydration_ftk(self):
        return self.temperature_problem.heat_of_hydration_ftk

    def get_E_alpha_fkt(self):
        return np.vectorize(self.mechanics_problem.E_fkt)

    def get_X_alpha_fkt(self):
        return self.mechanics_problem.general_hydration_fkt


class ConcreteTempHydrationModel(df.NonlinearProblem):
    def __init__(self, mesh, p, pv_name='temp_output', **kwargs):
        df.NonlinearProblem.__init__(self)  # apparently required to initialize things
        self.p = p

        if mesh != None:
            # initialize possible paraview output
            self.pv_file = df.XDMFFile(pv_name + '.xdmf')
            self.pv_file.parameters["flush_output"] = True
            self.pv_file.parameters["functions_share_mesh"] = True
            # function space for single value per element, required for plot of quadrature space values

            # initialize timestep, musst be reset using .set_timestep(dt)
            self.dt = 0
            self.dt_form = df.Constant(self.dt)

            if self.p.degree == 1:
                self.visu_space = df.FunctionSpace(mesh, "DG", 0)
            else:
                self.visu_space = df.FunctionSpace(mesh, "P", 1)

            metadata = {"quadrature_degree": self.p.degree, "quadrature_scheme": "default"}
            dxm = df.dx(metadata=metadata)

            # solution field
            self.V = df.FunctionSpace(mesh, 'P', self.p.degree)

            # generic quadrature function space
            cell = mesh.ufl_cell()
            q = "Quadrature"
            quadrature_element = df.FiniteElement(q, cell, degree=self.p.degree, quad_scheme="default")
            q_V = df.FunctionSpace(mesh, quadrature_element)

            # quadrature functions
            self.q_T = df.Function(q_V, name="temperature")
            self.q_alpha = df.Function(q_V, name="degree of hydration")
            self.q_alpha_n = df.Function(q_V, name="degree of hydration last time step")
            self.q_delta_alpha = df.Function(q_V, name="inrease in degree of hydration")
            self.q_ddalpha_dT = df.Function(q_V, name="derivative of delta alpha wrt temperature")

            # empfy list for newton iteration to compute delta alpha using the last value as starting point
            self.delta_alpha_n_list = np.full(np.shape(self.q_alpha_n.vector().get_local()), 0.2)
            # empfy list for newton iteration to compute delta alpha using the last value as starting point
            self.delta_alpha_guess = np.full(np.shape(self.q_alpha_n.vector().get_local()), 0.5)

            # scalars for the analysis of the heat of hydration
            self.alpha = 0
            self.delta_alpha = 0

            # Define variational problem
            self.T = df.Function(self.V)  # temperature
            self.T_n = df.Function(self.V)  # overwritten later...
            T_ = df.TrialFunction(self.V)  # temperature
            vT = df.TestFunction(self.V)

            # normal form
            R_ufl = df.Constant(self.p.vol_heat_cap) * (self.T) * vT * dxm
            R_ufl += self.dt_form * df.dot(df.Constant(self.p.themal_cond) * df.grad(self.T), df.grad(vT)) * dxm
            R_ufl += -  df.Constant(self.p.vol_heat_cap) * self.T_n * vT * dxm
            # quadrature point part

            self.R = R_ufl - df.Constant(
                self.p.Q_pot * self.p.density_binder * self.p.b_ratio) * self.q_delta_alpha * vT * dxm

            # derivative
            # normal form
            dR_ufl = df.derivative(R_ufl, self.T)
            # quadrature part
            self.dR = dR_ufl - df.Constant(
                self.p.Q_pot * self.p.density_binder * self.p.b_ratio) * self.q_ddalpha_dT * T_ * vT * dxm

            # setup projector to project continuous funtionspace to quadrature
            self.project_T = LocalProjector(self.T, q_V, dxm)

            self.assembler = None  # set as default, to check if bc have been added???

    def delta_alpha_fkt(self, delta_alpha, alpha_n, T):
        return delta_alpha - self.dt * self.affinity(delta_alpha, alpha_n) * self.temp_adjust(T)

    def delta_alpha_prime(self, delta_alpha, alpha_n, T):
        return 1 - self.dt * self.daffinity_ddalpha(delta_alpha, alpha_n) * self.temp_adjust(T)

    def heat_of_hydration_ftk(self, T, time_list, dt, parameter):

        def interpolate(x, x_list, y_list):
            # assuming ordered x list

            i = 0
            # check if x is in the dataset
            if x > x_list[-1]:
                print(' * Warning!!!: Extrapolation!!!')
                point1 = (x_list[-2], y_list[-2])
                point2 = (x_list[-1], y_list[-1])
            elif x < x_list[0]:
                print(' * Warning!!!: Extrapolation!!!')
                point1 = (x_list[0], y_list[0])
                point2 = (x_list[1], y_list[1])
            else:
                while x_list[i] < x:
                    i += 1
                point1 = (x_list[i - 1], y_list[i - 1])
                point2 = (x_list[i], y_list[i])

            slope = (point2[1] - point1[1]) / (point2[0] - point1[0])
            x_increment = x - point1[0]
            y_increment = slope * x_increment
            y = point1[1] + y_increment

            return y

        # get tmax, identify number of time steps, then interpolate data
        # assuming time list is ordered!!!
        tmax = time_list[-1]

        # set paramters
        self.p.B1 = parameter['B1']
        self.p.B2 = parameter['B2']
        self.p.eta = parameter['eta']
        self.p.alpha_max = parameter['alpha_max']
        self.p.E_act = parameter['E_act']
        self.p.T_ref = parameter['T_ref']
        self.p.Q_pot = parameter['Q_pot']

        # set time step
        self.dt = dt

        t = 0
        time = [0.0]
        heat = [0.0]
        alpha_list = [0.0]
        alpha = 0
        delta_alpha = 0.0

        error_flag = False
        while t < tmax:
            # compute delta_alpha
            try:
                delta_alpha = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha, T + self.p.zero_C),
                                                    fprime=self.delta_alpha_prime, x0=delta_alpha)
                if delta_alpha < 0:
                    raise Exception(
                        f'Problem with solving for delta alpha. Result is negative for starting delta alpha = {delta_alpha}')
            except:
                delta_alpha = 0.2
                try:
                    delta_alpha = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha, T + self.p.zero_C),
                                                        fprime=self.delta_alpha_prime, x0=delta_alpha)
                    if delta_alpha < 0:
                        raise Exception(
                            'Problem with solving for delta alpha. Result is negative for starting delta alpha = 0.2')
                except:
                    delta_alpha = 0.5
                    try:
                        delta_alpha = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha, T + self.p.zero_C),
                                                            fprime=self.delta_alpha_prime, x0=delta_alpha)
                        if delta_alpha < 0:
                            raise Exception(
                                'Problem with solving for delta alpha. Result is negative for starting delta alpha = 0.5')
                    except:
                        delta_alpha = 1.0

                        try:
                            delta_alpha = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha, T + self.p.zero_C),
                                                                fprime=self.delta_alpha_prime, x0=delta_alpha)
                            if delta_alpha < 0:
                                raise Exception('Problem with solving for delta alpha. Result is negative.')
                        except:
                            error_flag = True
                            break

            # update alpha
            alpha = delta_alpha + alpha
            # save heat of hydration
            alpha_list.append(alpha)
            heat.append(alpha * self.p.Q_pot)

            # timeupdate
            t = t + self.dt
            time.append(t)

        # if there was a probem with the computation (bad input values), return zero
        if error_flag:
            heat_interpolated = np.zeros_like(time_list)
            alpha_interpolated = np.zeros_like(time_list)
        else:
            # interpolate heat to match time_list
            heat_interpolated = []
            alpha_interpolated = []
            for value in time_list:
                heat_interpolated.append(interpolate(value, time, heat))
                alpha_interpolated.append(interpolate(value, time, alpha_list))

        return np.asarray(heat_interpolated) / 1000, np.asarray(alpha_interpolated)

    def get_affinity(self):
        alpha_list = []
        affinity_list = []
        for val in range(1000):
            alpha = val / 1000
            alpha_list.append(alpha)
            affinity_list.append(self.affinity(alpha, 0))

        return np.asarray(alpha_list), np.asarray(affinity_list)

    def evaluate_material(self):
        # project temperautre onto quadrature spaces
        self.project_T(self.q_T)

        # convert quadrature spaces to numpy vector
        temperature_list = self.q_T.vector().get_local()
        alpha_n_list = self.q_alpha_n.vector().get_local()

        # solve for alpha at each quadrature point
        # here the newton raphson method of the scipy package is used
        # the zero value of the delta_alpha_fkt is found for each entry in alpha_n_list is found. the corresponding temparature
        # is given in temperature_list and as starting point the value of last step used from delta_alpha_n
        try:
            delta_alpha_list = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha_n_list, temperature_list),
                                                     fprime=self.delta_alpha_prime, x0=self.delta_alpha_n_list)
            # I dont trust the algorithim!!! check if only applicable results are obtained
        except:
            # AAAAAAHHHH, negative delta alpha!!!!
            # NO PROBLEM!!!, different starting value!
            delta_alpha_list = scipy.optimize.newton(self.delta_alpha_fkt, args=(alpha_n_list, temperature_list),
                                                     fprime=self.delta_alpha_prime, x0=self.delta_alpha_guess)
            if np.any(delta_alpha_list < 0.0):
                print('AAAAAAHHHH, negative delta alpha!!!!')
                raise Exception(
                    'There is a problem with the alpha computation/initial guess, computed delta alpha is negative.')

        # save the delta alpha for next iteration as starting guess
        self.delta_alpha_n_list = delta_alpha_list

        # compute current alpha
        alpha_list = alpha_n_list + delta_alpha_list
        # compute derivative of delta alpha with respect to temperature for rhs
        ddalpha_dT_list = self.dt * self.affinity(alpha_list, alpha_n_list) * self.temp_adjust_tangent(temperature_list)

        # project lists onto quadrature spaces
        set_q(self.q_alpha, alpha_list)
        set_q(self.q_delta_alpha, delta_alpha_list)
        set_q(self.q_ddalpha_dT, ddalpha_dT_list)

    def update_history(self):
        self.T_n.assign(self.T)  # save temparature field
        self.q_alpha_n.assign(self.q_alpha)  # save alpha field

    def set_timestep(self, dt):
        self.dt = dt
        self.dt_form.assign(df.Constant(self.dt))

    def set_initial_T(self, T):
        # set initial temperature, in kelvin
        T0 = df.Expression('t_zero', t_zero=T + self.p.zero_C, degree=0)
        self.T_n.interpolate(T0)
        self.T.interpolate(T0)

    def set_bcs(self, bcs):
        # Only now (with the bcs) can we initialize the assembler
        self.assembler = df.SystemAssembler(self.dR, self.R, bcs)

    def F(self, b, x):
        if self.dt <= 0:
            raise RuntimeError("You need to `.set_timestep(dt)` larger than zero before the solve!")
        if not self.assembler:
            raise RuntimeError("You need to `.set_bcs(bcs)` before the solve!")
        self.evaluate_material()
        self.assembler.assemble(b, x)

    def J(self, A, x):
        self.assembler.assemble(A)

    def pv_plot(self, t=0):
        # paraview export

        # temperature plot
        T_plot = df.project(self.T, self.V)
        T_plot.rename("Temperature", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(T_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

        # degree of hydration plot
        alpha_plot = df.project(self.q_alpha, self.visu_space,
                                form_compiler_parameters={'quadrature_degree': self.p.degree})
        alpha_plot.rename("DOH", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(alpha_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

    def temp_adjust(self, T):
        val = 1
        if self.p.temp_adjust_law == 'exponential':
            val = np.exp(-self.p.E_act / self.p.igc * (1 / T - 1 / (self.p.T_ref + self.p.zero_C)))
        elif self.p.temp_adjust_law == 'off':
            pass
        else:
            # TODO throw correct error
            raise Exception(
                f'Warning: Incorrect temp_adjust_law {self.p.temp_adjust_law} given, only "exponential" and "off" implemented')
        return val

        # derivative of the temperature adjustment factor with respect to the temperature

    def temp_adjust_tangent(self, T):
        val = 0
        if self.p.temp_adjust_law == 'exponential':
            val = self.p.E_act / self.p.igc / T ** 2
        return val

    # affinity function
    def affinity(self, delta_alpha, alpha_n):
        affinity = self.p.B1 * (self.p.B2 / self.p.alpha_max + delta_alpha + alpha_n) * (
                self.p.alpha_max - (delta_alpha + alpha_n)) * np.exp(
            -self.p.eta * (delta_alpha + alpha_n) / self.p.alpha_max)
        return affinity

    # derivative of affinity with respect to delta alpha
    def daffinity_ddalpha(self, delta_alpha, alpha_n):
        affinity_prime = self.p.B1 * np.exp(-self.p.eta * (delta_alpha + alpha_n) / self.p.alpha_max) * (
                (self.p.alpha_max - (delta_alpha + alpha_n)) * (
                self.p.B2 / self.p.alpha_max + (delta_alpha + alpha_n)) * (
                        -self.p.eta / self.p.alpha_max) - self.p.B2 / self.p.alpha_max - 2 * (
                        delta_alpha + alpha_n) + self.p.alpha_max)
        return affinity_prime


class ConcreteMechanicsModel(df.NonlinearProblem):
    def __init__(self, mesh, p, pv_name='mechanics_output', **kwargs):
        df.NonlinearProblem.__init__(self)  # apparently required to initialize things
        self.p = p

        if self.p.dim == 1:
            self.stress_vector_dim = 1
        elif self.p.dim == 2:
            self.stress_vector_dim = 3
        elif self.p.dim == 3:
            self.stress_vector_dim = 6

        # todo: I do not like the "meshless" setup right now
        if mesh != None:
            # initialize possible paraview output
            self.pv_file = df.XDMFFile(pv_name + '.xdmf')
            self.pv_file.parameters["flush_output"] = True
            self.pv_file.parameters["functions_share_mesh"] = True
            # function space for single value per element, required for plot of quadrature space values

            #
            if self.p.degree == 1:
                self.visu_space = df.FunctionSpace(mesh, "DG", 0)
                self.visu_space_T = df.TensorFunctionSpace(mesh, "DG", 0)
            else:
                self.visu_space = df.FunctionSpace(mesh, "P", 1)
                self.visu_space_T = df.TensorFunctionSpace(mesh, "P", 1)

            metadata = {"quadrature_degree": self.p.degree, "quadrature_scheme": "default"}
            dxm = df.dx(metadata=metadata)

            # solution field
            self.V = df.VectorFunctionSpace(mesh, 'P', self.p.degree)

            # generic quadrature function space
            cell = mesh.ufl_cell()
            q = "Quadrature"

            quadrature_element = df.FiniteElement(q, cell, degree=self.p.degree, quad_scheme="default")
            quadrature_vector_element = df.VectorElement(q, cell, degree=self.p.degree, dim=self.stress_vector_dim,
                                                         quad_scheme="default")
            q_V = df.FunctionSpace(mesh, quadrature_element)
            q_VT = df.FunctionSpace(mesh, quadrature_vector_element)

            # quadrature functions
            self.q_E = df.Function(q_V, name="youngs modulus")
            self.q_fc = df.Function(q_V, name="compressive strength")
            self.q_ft = df.Function(q_V, name="tensile strength")
            self.q_yield = df.Function(q_V, name="yield criterion")
            self.q_alpha = df.Function(q_V, name="degree of hydration")

            self.q_sigma = df.Function(q_VT, name="stress")
            # initialize degree of hydration to 1, in case machanics module is run without hydration coupling
            self.q_alpha.vector()[:] = 1

            # Define variational problem
            self.u = df.Function(self.V)  # displacement
            v = df.TestFunction(self.V)

            # Elasticity parameters without multiplication with E
            x_mu = 1.0 / (2.0 * (1.0 + self.p.nu))
            x_lambda = 1.0 * self.p.nu / ((1.0 + self.p.nu) * (1.0 - 2.0 * self.p.nu))

            # Stress computation for linear elastic problem without multiplication with E
            def x_sigma(v):
                return 2.0 * x_mu * df.sym(df.grad(v)) + x_lambda * df.tr(df.sym(df.grad(v))) * df.Identity(len(v))

            # Volume force            
            if self.p.dim == 1:
                f = df.Constant(-self.p.g * self.p.density)
            elif self.p.dim == 2:
                f = df.Constant((0, -self.p.g * self.p.density))
            elif self.p.dim == 3:
                f = df.Constant((0, 0, -self.p.g * self.p.density))

            self.sigma_ufl = self.q_E * x_sigma(self.u)

            R_ufl = self.q_E * df.inner(x_sigma(self.u), df.sym(df.grad(v))) * dxm
            R_ufl += - df.inner(f, v) * dxm  # add volumetric force, aka gravity (in this case)
            # quadrature point part
            self.R = R_ufl  # - Constant(p.Q_inf) * self.q_delta_alpha * vT * dxm

            # derivative
            # normal form
            dR_ufl = df.derivative(R_ufl, self.u)
            # quadrature part
            self.dR = dR_ufl  # - Constant(p.Q_inf) * self.q_ddalpha_dT * T_ * vT * dxm

            self.project_sigma = LocalProjector(self.sigma_voigt(self.sigma_ufl), q_VT, dxm)

            self.assembler = None  # set as default, to check if bc have been added???

    def sigma_voigt(self, s):
        # 1D option
        if s.ufl_shape == (1, 1):
            stress_vector = df.as_vector((s[0, 0]))
        # 2D option
        elif s.ufl_shape == (2, 2):
            stress_vector = df.as_vector((s[0, 0], s[1, 1], s[0, 1]))
        # 3D option
        elif s.ufl_shape == (3, 3):
            stress_vector = df.as_vector((s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[0, 2]))
        else:
            raise ('Problem with stress tensor shape for voigt notation')

        return stress_vector

    def E_fkt(self, alpha, parameters):

        if alpha < parameters['alpha_t']:
            E = parameters['E_inf'] * alpha / parameters['alpha_t'] * (
                        (parameters['alpha_t'] - parameters['alpha_0']) / (1 - parameters['alpha_0'])) ** parameters[
                    'a_E']
        else:
            E = parameters['E_inf'] * ((alpha - parameters['alpha_0']) / (1 - parameters['alpha_0'])) ** parameters[
                'a_E']
        return E

    def general_hydration_fkt(self, alpha, parameters):

        return parameters['X_inf'] * alpha ** (parameters['a_X'])

    def principal_stress(self, stresses):
        # checking type of problem
        n = stresses.shape[1]  # number of stress components in stress vector
        # finding eigenvalues of symmetric stress tensor
        # 1D problem
        if n == 1:
            principal_stresses = stresses
        # 2D problem
        elif n == 3:
            # the following uses
            # lambda**2 - tr(sigma)lambda + det(sigma) = 0, solve for lambda using pq formula
            p = - (stresses[:, 0] + stresses[:, 1])
            q = stresses[:, 0] * stresses[:, 1] - stresses[:, 2] ** 2

            D = p ** 2 / 4 - q  # help varibale
            assert np.all(D >= -1.0e-15)  # otherwise problem with imaginary numbers
            sqrtD = np.sqrt(D)

            eigenvalues_1 = -p / 2.0 + sqrtD
            eigenvalues_2 = -p / 2.0 - sqrtD

            # strack lists as array
            principal_stresses = np.column_stack((eigenvalues_1, eigenvalues_2))

            # principal_stress = np.array([ev1p,ev2p])
        elif n == 6:
            # for a symetric stress vector a b c e f d we need to solve:
            # x**3 - x**2(a+b+c) - x(e**2+f**2+d**2-ab-bc-ac) + (abc-ae**2-bf**2-cd**2+2def) = 0, solve for x
            principal_stresses = np.empty([len(stresses), 3])
            # currently slow solution with loop over all stresses and subsequent numpy function call:
            for i, stress in enumerate(stresses):
                # convert voigt to tensor, (00,11,22,12,02,01)
                stress_tensor = np.zeros((3, 3))
                stress_tensor[0][0] = stress[0]
                stress_tensor[1][1] = stress[1]
                stress_tensor[2][2] = stress[2]
                stress_tensor[0][1] = stress[5]
                stress_tensor[1][2] = stress[3]
                stress_tensor[0][2] = stress[4]
                stress_tensor[1][0] = stress[5]
                stress_tensor[2][1] = stress[3]
                stress_tensor[2][0] = stress[4]
                # use numpy for eigenvalues
                principal_stress = np.linalg.eigvalsh(stress_tensor)
                # sort principal stress from lagest to smallest!!!
                principal_stresses[i] = -np.sort(-principal_stress)

        return principal_stresses

    def yield_surface(self, stresses, ft, fc):
        # function for approximated yield surface
        # first approximation, could be changed if we have numbers/information
        fc2 = fc
        # pass voigt notation and compute the principal stress
        p_stresses = self.principal_stress(stresses)

        # get the principle tensile stresses
        t_stresses = np.where(p_stresses < 0, 0, p_stresses)

        # get dimension of problem, ie. length of list with principal stresses
        n = p_stresses.shape[1]
        # check case
        if n == 1:
            # rankine for the tensile region
            rk_yield_vals = t_stresses[:, 0] - ft[:]

            # invariants for drucker prager yield surface
            I1 = stresses[:, 0]
            I2 = np.zeros_like(I1)
        # 2D problem
        elif n == 2:

            # rankine for the tensile region
            rk_yield_vals = (t_stresses[:, 0] ** 2 + t_stresses[:, 1] ** 2) ** 0.5 - ft[:]

            # invariants for drucker prager yield surface
            I1 = stresses[:, 0] + stresses[:, 1]
            I2 = ((stresses[:, 0] + stresses[:, 1]) ** 2 - ((stresses[:, 0]) ** 2 + (stresses[:, 1]) ** 2)) / 2

        # 3D problem
        elif n == 3:
            # rankine for the tensile region
            rk_yield_vals = (t_stresses[:, 0] ** 2 + t_stresses[:, 1] ** 2 + t_stresses[:, 2] ** 2) ** 0.5 - ft[:]

            # invariants for drucker prager yield surface
            I1 = stresses[:, 0] + stresses[:, 1] + stresses[:, 2]
            I2 = ((stresses[:, 0] + stresses[:, 1] + stresses[:, 2]) ** 2 - (
                        (stresses[:, 0]) ** 2 + (stresses[:, 1]) ** 2 + (stresses[:, 2]) ** 2)) / 2
        else:
            raise ('Problem with input to yield surface, the array with stress values has the wrong size ')

        J2 = 1 / 3 * I1 ** 2 - I2
        beta = (3.0 ** 0.5) * (fc2 - fc) / (2 * fc2 - fc)
        Hp = fc2 * fc / ((3.0 ** 0.5) * (2 * fc2 - fc))

        dp_yield_vals = beta / 3 * I1 + J2 ** 0.5 - Hp

        # TODO: is this "correct", does this make sense? for a compression state, what if rk yield > dp yield???
        yield_vals = np.maximum(rk_yield_vals, dp_yield_vals)

        return np.asarray(yield_vals)

    def evaluate_material(self):
        # convert quadrature spaces to numpy vector
        alpha_list = self.q_alpha.vector().get_local()

        parameters = {}
        parameters['alpha_t'] = self.p.alpha_t
        parameters['E_inf'] = self.p.E_28
        parameters['alpha_0'] = self.p.alpha_0
        parameters['a_E'] = self.p.a_E
        # vectorize the function for speed up
        E_fkt_vectorized = np.vectorize(self.E_fkt)
        E_list = E_fkt_vectorized(alpha_list, parameters)

        parameters = {}
        parameters['X_inf'] = self.p.fc_inf
        parameters['a_X'] = self.p.a_fc

        fc_list = self.general_hydration_fkt(alpha_list, parameters)

        parameters = {}
        parameters['X_inf'] = self.p.ft_inf
        parameters['a_X'] = self.p.a_ft

        ft_list = self.general_hydration_fkt(alpha_list, parameters)

        # now do the yield function thing!!!
        # I need stresses!!!
        # get stress values
        self.project_sigma(self.q_sigma)

        sigma_list = self.q_sigma.vector().get_local().reshape((-1, self.stress_vector_dim))

        # compute the yield values (values > 0 : failure)
        yield_list = self.yield_surface(sigma_list, ft_list, fc_list)

        # # project lists onto quadrature spaces
        set_q(self.q_E, E_list)
        set_q(self.q_fc, fc_list)
        set_q(self.q_ft, ft_list)
        set_q(self.q_yield, yield_list)

    def update_history(self):
        # no history field currently
        pass

    def set_timestep(self, dt):
        self.dt = dt
        self.dt_form.assign(df.Constant(self.dt))

    def set_bcs(self, bcs):
        # Only now (with the bcs) can we initialize the assembler
        self.assembler = df.SystemAssembler(self.dR, self.R, bcs)

    def F(self, b, x):
        # if self.dt <= 0:
        #    raise RuntimeError("You need to `.set_timestep(dt)` larger than zero before the solve!")
        if not self.assembler:
            raise RuntimeError("You need to `.set_bcs(bcs)` before the solve!")
        self.evaluate_material()
        self.assembler.assemble(b, x)

    def J(self, A, x):
        self.assembler.assemble(A)

    def pv_plot(self, t=0):
        # paraview export

        # displacement plot
        u_plot = df.project(self.u, self.V)
        u_plot.rename("Displacement", "test string, what does this do??")  # TODO: what does the second string do?
        self.pv_file.write(u_plot, t, encoding=df.XDMFFile.Encoding.ASCII)

        # Elasticity parameters without multiplication with E
        x_mu = 1.0 / (2.0 * (1.0 + self.p.nu))
        x_lambda = 1.0 * self.p.nu / ((1.0 + self.p.nu) * (1.0 - 2.0 * self.p.nu))

        def x_sigma(v):
            return 2.0 * x_mu * df.sym(df.grad(v)) + x_lambda * df.tr(df.sym(df.grad(v))) * df.Identity(len(v))

        sigma_plot = df.project(self.sigma_ufl, self.visu_space_T,
                                form_compiler_parameters={'quadrature_degree': self.p.degree})
        E_plot = df.project(self.q_E, self.visu_space, form_compiler_parameters={'quadrature_degree': self.p.degree})
        fc_plot = df.project(self.q_fc, self.visu_space, form_compiler_parameters={'quadrature_degree': self.p.degree})
        ft_plot = df.project(self.q_ft, self.visu_space, form_compiler_parameters={'quadrature_degree': self.p.degree})
        yield_plot = df.project(self.q_yield, self.visu_space,
                                form_compiler_parameters={'quadrature_degree': self.p.degree})
        #
        E_plot.rename("Young's Modulus", "test string, what does this do??")  # TODO: what does the second string do?
        fc_plot.rename("Compressive strength",
                       "test string, what does this do??")  # TODO: what does the second string do?
        ft_plot.rename("Tensile strength", "test string, what does this do??")  # TODO: what does the second string do?
        yield_plot.rename("Yield surface", "test string, what does this do??")  # TODO: what does the second string do?
        sigma_plot.rename("Stress", "test string, what does this do??")  # TODO: what does the second string do?

        self.pv_file.write(E_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
        self.pv_file.write(fc_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
        self.pv_file.write(ft_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
        self.pv_file.write(yield_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
        self.pv_file.write(sigma_plot, t, encoding=df.XDMFFile.Encoding.ASCII)
