import fenics_concrete
import numpy as np
from abc import ABC, abstractmethod

# TODO: inherit from abstract baselin class for solvers

class ForwardBase(ABC):
    """Base class for forward solvers
    """
    def __init__(self):
        pass

    @abstractmethod
    def solve(self,latents:list,inp_solver:dict, **kwargs)->list:
        pass



class HydrationSolverWrapper(ForwardBase):
    def __init__(self):
        super().__init__()
    def solve(self,latents:list,inp_solver:dict, **kwargs)->list:
        parameter = fenics_concrete.Parameters()  # using the current default values
        # -- latents -----
        # parameter['B1'] = 2.916E-4  # in 1/s (le 0, < 0.1)
        # parameter['B2'] = 0.0024229  # - (le 0, smaller 1)
        # parameter['eta'] = 5.554  # something about diffusion (should be larger 0)
        # parameter['T_ref'] = 25  # reference temperature in degree celsius
        # parameter['Q_pot'] = 500e3 # potential heat per weight of binder in J/kg

        # -- adding scaling back the values
        parameter['B1'] = latents[0]*1e-04  # in 1/s (le 0, < 0.1)
        parameter['B2'] = latents[1]*1e-03  # - (le 0, smaller 1)
        parameter['eta'] = latents[2]  # something about diffusion (should be larger 0)
        parameter['Q_pot'] = latents[3]*1e05  # potential heat per weight of binder in J/kg

        # -- observed inputs
        parameter['igc'] = 8.3145  # ideal gas constant in [J/K/mol], CONSTANT!!!
        parameter['zero_C'] = 273.15  # in Kelvin, CONSTANT!!!
        parameter['E_act'] = 47002  # activation energy in Jmol^-1 (no relevant limits) (Depends only on simulated temp, if that is not change no need to infer E_act)
        parameter['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c (larger 0 and max 1)
        parameter['T_ref'] = 25  # reference temperature in degree celsius

        # this is the minimal time step used in the simulation
        # using a larger value will increase the speed but decrease the accuracy
        dt = 300 # value in seconds

        # this is the simulated temperature, needs to be adjusted depending on the temperature of the experimental data
        T = inp_solver['T_rxn'] # can be 20,40,60 as pert the exp values
        # this is the list of measured time data as given by the experiments
        #time_list = [0,5000,10000,20000,100000]
        time_list = inp_solver['time_list']

        # initiate material problem, for this the "fenics_concrete" conda package needs to be installed
        # use: 'mamba install -c etamsen fenics_concrete"
        problem = fenics_concrete.ConcreteThermoMechanical()

        # get the hydration function
        # this might change in the future to make it more easily accessible but for now it should work like this
        hydration_fkt = problem.get_heat_of_hydration_ftk()
        # the results are a heat list and a degree of hydration list, which you can ignore for now
        heat_list, doh_list= hydration_fkt(T, time_list, dt, parameter)

        return heat_list

# Homogenization solver
class HomogenizationSolverWrapper(ForwardBase):
    def __init__(self):
        super().__init__()
    def solve(self,latents:list,inp_solver:dict, **kwargs)->list:
        return NotImplementedError

# write pytests
def test_hydration_solver_wrapper():
    # -- observed inputs
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = [0,5000,10000,20000,100000]

    # -- latents -----
    b = [2.916,2.4229,5.554,5]
    hydration_solver = HydrationSolverWrapper()
    heat_list = hydration_solver.solve(latents=b,inp_solver=inp_solver)
    #heat_list = hydration_solver_wrapper(b,inp_solver)
    print(f'heat_list = {heat_list}')

    # -- expected outputs
    heat_list_exp =[  0.        ,   2.14549938,   7.1823244 ,  34.34254352,
       233.33527714]
    # assert the values are approximately equal
    # write assert statement also
    assert np.allclose(heat_list,heat_list_exp,atol=1e-3), "The heat list is not equal to the expected values"

if __name__ == "__main__":
    test_hydration_solver_wrapper()

