import numpy as np
import os

from lebedigital.demonstrator_calibration.prior import prior
from lebedigital.unit_registry import ureg
import torch as th
th.set_default_dtype(th.float64)

def computation_paste_strength_stiffness(slag_ratio:float,gaussian_mean:str, gaussian_cov:str, seed:int):
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
    paste_youngs_modulus : float / pint stress unit, will be in 'GPa'
        approximated youngs modulus of paste
    paste_strength : float / pint stress unit, will be in 'MPa'
        approximated compressive strength of paste
    """

    assert os.path.isfile(gaussian_mean), f"The file {gaussian_mean} does not exist"

    mean_model = th.jit.load(gaussian_mean)

    # load cov parameters
    # ----- check if the file exists
    assert os.path.isfile(gaussian_cov), f"The file {gaussian_cov} does not exist"
    cov = np.genfromtxt(gaussian_cov, delimiter=',')
    assert len(cov) == 3, f"The cov matrix should be 2x2, but is is not"

    

    # define the distribution
    # ----- set the seed
    th.manual_seed(seed=seed)
    # ----- define the transformation, ugly now, should be imported from somewhere
    def transformed_back(samples):
        shape = samples.shape
        samples[:, 0] = samples[:, 0] * 1e09
        samples[:, 1] = samples[:, 1]* 1e07
        assert samples.shape == shape, "shape of the samples is changed"
        return samples.flatten()
    
    # ----- sample from the distribution
    dist = prior(mean=mean_model, cov_params=cov, cov_type='full',latent_dim=2)
    sample = dist.sample(x=[slag_ratio],n_samples=1)
    E_paste, fc_paste = transformed_back(sample) # returns in Pa

    # converting to Gpa and Mpa respectively
    E_paste = E_paste*1e-09
    fc_paste = fc_paste*1e-06

    return E_paste * ureg('GPa') , fc_paste * ureg('MPa')