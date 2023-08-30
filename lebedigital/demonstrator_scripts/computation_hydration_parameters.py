import numpy as np
import os
from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_calibration.prior import prior
import torch as th
import numpy as np

th.set_default_dtype(th.float64)

def computation_hydration_parameters(slag_ratio:float,gaussian_mean:str, gaussian_cov:str, seed:int):
    """_summary_

    Parameters
    ----------
    slag_ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    gaussian_mean : str
        location to the .pt file of the pretrained NN. This should be pre-scripted torch model, meaning 
        need to define the class. 
    gaussian_cov : str
        location to the .csv file with the cov.
    seed : int
        The value of the seed which is passed to the prior function.
    
    Returns
    -------
    B1 : float / pint unit will be in '1/s'
        hydration parameter
    B2 : float / pint unitless unit
        hydration parameter
    eta : float / pint unitless unit
        hydration parameter
    E_act : float / pint unit, will be in 'J/mol'
        activation energy - hydration parameter
    T_ref : float / pint unit, will be in 'degree Celsius'
        reference temperature - hydration parameter
    Q_pot : float / pint unit, will be in 'J/kg'
        maximum potential hydration parameter
    """

    # load the Neural Network
    # ----- check if the file exists
    assert os.path.isfile(gaussian_mean), f"The file {gaussian_mean} does not exist"

    mean_model = th.jit.load(gaussian_mean)

    # load cov parameters
    # ----- check if the file exists
    assert os.path.isfile(gaussian_cov), f"The file {gaussian_cov} does not exist"
    cov = np.genfromtxt(gaussian_cov, delimiter=',')
    assert len(cov) == 10, f"The cov matrix should be 4x4, but is is not"

    # define the distribution
    # ----- set the seed
    th.manual_seed(seed=seed)
    # ----- define the transformation, ugly now, should be imported for somewhere
    def transformed_back(samples):
        shape = samples.shape
        # exp transform to the last dimention
        samples[:,0] = samples[:,0] * 1e-04
        samples[:,1] = np.exp(samples[:,1]) # the 3rd value (eta) is not scaled
        samples[:,3] = samples[:,3] * 1e05
        assert samples.shape == shape, "shape of the samples is changed"
        return samples.flatten()
    
    # ----- sample from the distribution
    dist = prior(mean=mean_model, cov_params=cov, cov_type='full',latent_dim=4)
    sample = dist.sample(x=[slag_ratio],n_samples=1)
    B1,B2,eta,Q_pot = transformed_back(sample)


    # other fixed parameters
    E_act = 5653 * 8.3145 * ureg('J/mol')  # activation energy in Jmol^-1
    Q_ = ureg.Quantity
    T_ref = Q_(20, ureg.degC) # this is 20 as the model-learning was done for 20 degC. This renders E_act useless in a way


    return B1 * ureg('1/s'), B2 * ureg(''), eta * ureg(''), \
           E_act , Q_pot * ureg('J/kg'), Q_(T_ref,ureg.degC)

if __name__ == "__main__":
    NN_path = 'usecases/optimization_paper/optimization_workflow/Inputs/NN_model_hydration_final.pt'
    cov_path = 'usecases/optimization_paper/optimization_workflow/Inputs/cov_parameters_hydration_final.csv'
    print(computation_hydration_parameters(0.2,NN_path,cov_path,seed=5))