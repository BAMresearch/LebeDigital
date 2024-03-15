import fenics_concrete
import numpy as np
from abc import ABC, abstractmethod
from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg

# TODO: inherit from abstract baselin class for solvers

class ForwardBase(ABC):
    """Base class for forward solvers
    """
    def __init__(self):
        pass
    @abstractmethod
    def _scale_back(self,latent:float):
        pass
    @abstractmethod
    def solve(self,latents:list,inp_solver:dict, **kwargs)->list:
        pass



class HydrationSolverWrapper(ForwardBase):
    def __init__(self):
        super().__init__()
    def _scale_back(self,latent:list):
        """this assumes inputs are log scaled, so here they are scaled back"""
        #return np.exp(latent)
        # expectes normalized values, so scale back to original values
        #std = np.array([1.9956, 247.6045,   1.8181,   2.5245]) # values from scipy optimizer
        #mean = np.array([  2.8128, 124.1033,   3.4967,   3.6444]) 
        # 'B_1', 'B_2', 'eta', 'Q_pot' assumes this order
        latent[1] = np.exp(latent[1])
        latent[-1] = np.exp(latent[-1]) # as E_a is always positive
        #latent_scaled_back = np.array(latent)*std + mean
        latent_scaled_back = np.array(latent)
        return latent_scaled_back
    def solve(self,latents:list,inp_solver:dict, **kwargs)->list:
        parameter = fenics_concrete.Parameters()  # using the current default values
        # -- latents -----
        # parameter['B1'] = 2.916E-4  # in 1/s (le 0, < 0.1)
        # parameter['B2'] = 0.0024229  # - (le 0, smaller 1)
        # parameter['eta'] = 5.554  # something about diffusion (should be larger 0)
        # parameter['T_ref'] = 25  # reference temperature in degree celsius
        # parameter['Q_pot'] = 500e3 # potential heat per weight of binder in J/kg

        # -- adding scaling back the values
        latent_scaled_back = self._scale_back(latents)
        parameter['B1'] = latent_scaled_back[0]*1e-04  # in 1/s (le 0, < 0.1)
        #parameter['B2'] = latent_scaled_back[1]*1e-03  # - (le 0, smaller 1)
        parameter['B2'] = latent_scaled_back[1]
        parameter['eta'] = latent_scaled_back[2]  # something about diffusion (should be larger 0)
        parameter['Q_pot'] = latent_scaled_back[3]*1e05  # potential heat per weight of binder in J/kg
        parameter['E_act'] = latent_scaled_back[4]*1e04  # activation energy in Jmol^-1 (no relevant limits) (Depends only on simulated temp, if that is not change no need to infer E_act)

        # -- scaling back the values
        # parameter['B1'] = self._scale_back(latents[0]) # in 1/s (le 0, < 0.1)
        # parameter['B2'] = self._scale_back(latents[1])  # - (le 0, smaller 1)
        # parameter['eta'] = self._scale_back(latents[2])  # something about diffusion (should be larger 0)
        # parameter['Q_pot'] = self._scale_back(latents[3])  # potential heat per weight of binder in J/kg

        # -- observed inputs
        parameter['igc'] = 8.3145  # ideal gas constant in [J/K/mol], CONSTANT!!!
        parameter['zero_C'] = 273.15  # in Kelvin, CONSTANT!!!
        #parameter['E_act'] = 47002  # activation energy in Jmol^-1 (no relevant limits) (Depends only on simulated temp, if that is not change no need to infer E_act)
        parameter['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c (larger 0 and max 1)
        #parameter['T_ref'] = inp_solver['T_rxn']  # reference temperature in degree celsius, if its = T_rxn, then E_ect doesnt matter
        #parameter['T_ref'] = 20
        parameter['T_ref'] = 22 # the temp the model learning was done on, this needs to bethe same later on also. 

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
    def _scale_back(self,latent:list):
        return [latent[0]*1e09,latent[1]*1e07]
    def solve(self,latents:list,inp_solver:dict=None, **kwargs)->list:
        # initialize dictionary
        parameters = {}

        # values taken from input/materials.json
        # paste data
        latent_scaled_back = self._scale_back(latents)
        parameters['paste_E'] = latent_scaled_back[0]* ureg('Pa')     #30e9 * ureg('Pa')
        parameters['paste_fc'] = latent_scaled_back[1]*ureg('Pa')     #30e6 * ureg('Pa')

        parameters['paste_nu'] = 0.26 * ureg('dimensionless')
        parameters['paste_C'] = 800 * ureg('J/kg/K')  # Specific Heat Capacity
        parameters['paste_kappa'] = 1.15 * ureg('W/m/K')  # Thermal conductivity
        parameters['paste_rho'] = 3100 * ureg('kg/m^3')
        parameters['paste_Q'] = 250000 * ureg('J/kg')

        # aggregate data
        parameters['aggregates_E'] = 65e9 * ureg('Pa')
        parameters['aggregates_nu'] = 0.25* ureg('dimensionless')
        parameters['aggregates_C'] = 800 * ureg('J/kg/K')  # Specific Heat Capacity
        parameters['aggregates_kappa'] = 3.1 * ureg('W/m/K')  # Thermal conductivity
        parameters['aggregates_rho'] = 2700 * ureg('kg/m^3')
        parameters['aggregates_vol_frac'] = 0.7* ureg('dimensionless')

        results = concrete_homogenization(parameters)
        result_for_learning = [results['E'].magnitude, results['fc'].magnitude] # both are in Pa
        return np.array(result_for_learning)






# write pytests
def test_hydration_solver_wrapper():
    # -- observed inputs
    inp_solver = {}
    inp_solver['T_rxn'] = 20
    inp_solver['time_list'] = [0,5000,10000,20000,100000]

    # -- latents -----
    # b = np.array([2.916,2.4229,5.554,5])
    # std = np.array([1.9956, 247.6045,   1.8181,   2.5245]) 
    # mean = np.array([  2.8128, 124.1033,   3.4967,   3.6444])
    # b = (b-mean)/std
    
    b = np.array([  2.91, np.log(2.422e-03),   3.4967,   3.6444, 4.7002])


    hydration_solver = HydrationSolverWrapper()
    heat_list = hydration_solver.solve(latents=b,inp_solver=inp_solver)
    #heat_list = hydration_solver_wrapper(b,inp_solver)
    print(f'heat_list = {heat_list}')

    # -- expected outputs
    heat_list_exp =[  0.,        17.61763829,  84.5571727, 181.80505507, 301.89535938]
    # assert the values are approximately equal
    # write assert statement also
    assert np.allclose(heat_list,heat_list_exp,atol=1e-3), "The heat list is not equal to the expected values"

def test_homogenization_solver():
    latents = [30,3]
    homogenization_solver = HomogenizationSolverWrapper()
    result = homogenization_solver.solve(latents=latents)
    print(f'result = {result}')
    result_correct = [51082128028.566986, 38101522.84263957]
    assert np.allclose(result,result_correct,atol=1e-3), "The homogenization solver is not working properly"

if __name__ == "__main__":
    test_hydration_solver_wrapper()
    #test_homogenization_solver()

