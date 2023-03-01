import numpy as np
from lebedigital.unit_registry import ureg
import torch as th
import numpy as np

th.set_default_dtype(th.float64)


@ureg.check('', '', '', '')
def dummy_hydration_parameters(slag_ratio, phi_mean: list, phi_cov: list, seed: int):
    """
    This is a dummy function to make the snakemake workflow work until the real function is ready
    It changes arbitrarily chosen values depending on slag content, not based on physics or anything.

    Parameters
    ----------
    seed : int
    phi_cov : list
    phi_mean : list
    slag_ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    phi_hydration: ??
        input from Atuls parameter identification

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

    # paste_youngs_modulus_min = 15 * ureg('GPa')
    # paste_youngs_modulus_max = 25 * ureg('GPa')
    # paste_youngs_modulus = paste_youngs_modulus_min + (paste_youngs_modulus_max-paste_youngs_modulus_min)*slag_ratio
    #
    # paste_strength_min = 8 * ureg('MPa')
    # paste_strength_max = 13 * ureg('MPa')
    # paste_strength = paste_strength_min + (paste_strength_max-paste_strength_min)*slag_ratio

    B1 = 2.916E-4 * ureg('1/s')  # in 1/s
    B2 = 0.0024229 * ureg('')  # -
    eta = 5.554 * ureg('')  # something about diffusion
    E_act = 5653 * 8.3145 * ureg('J/mol')  # activation energy in Jmol^-1
    Q_ = ureg.Quantity
    T_ref = Q_(25, ureg.degC)

    Q_pot_min = 100000 * ureg('J/kg')
    Q_pot_max = 300000 * ureg('J/kg')
    Q_pot = Q_pot_max - (Q_pot_max - Q_pot_min) * slag_ratio

    # ATUL : temporary q(b|x)~N(mu,cov), b=(B1,B2,eta,E_ect,Q_pot,T_ref), x = slag ratio
    th.manual_seed(seed=seed)
    no_of_parameter = 6
    assert len(phi_mean) == no_of_parameter
    assert np.array(phi_mean).ndim == 2

    mean = np.matmul(np.array(phi_mean)[:, :-1], np.atleast_1d(slag_ratio)) + np.array(phi_mean)[:,
                                                                              -1]  # assuming linear
    dist = th.distributions.MultivariateNormal(loc=th.as_tensor(mean), covariance_matrix=th.tensor(phi_cov))
    B1, B2, eta, E_act, Q_pot, T_ref = dist.sample()

    return B1.item() * ureg('1/s'), B2.item() * ureg(''), eta.item() * ureg(''), \
           E_act.item() * ureg('J/mol'), Q_pot.item() * ureg('J/kg'), Q_(T_ref.item(),ureg.degC)


if __name__ == "__main__":
    # test while developing this
    # slight slope for for all except Q. Q follows the same relation as Erik had proposed
    # Q_pot = Q_pot_max - (Q_pot_max - Q_pot_min) * slag_ratio
    # also guessing 5% noise. sigma^2 = (5/100)^2 = 0.0025

    phi_mean = [[0.01 * 2.916E-4, 2.916E-4], [0.01 * 0.0024229, 0.0024229], [0.01 * 5.554, 5.554],
                [0.01 * 5653 * 8.3145, 5653 * 8.3145], [-200000, 300000], [0.01 * 25, 25]]
    phi_cov = np.diag(0.0025 * np.array(phi_mean)[:, 1]).tolist()
    seed = 7
    print(dummy_hydration_parameters(0.5, phi_mean=phi_mean, phi_cov=phi_cov, seed=seed))
