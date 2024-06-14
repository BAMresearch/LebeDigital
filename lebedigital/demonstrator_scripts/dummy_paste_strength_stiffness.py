import numpy as np

from lebedigital.unit_registry import ureg
import torch as th
th.set_default_dtype(th.float64)


@ureg.check('', '', '', '')
def dummy_paste_strength_stiffness(slag_ratio, phi_mean, phi_cov, seed):
    """
    This is a dummy function to make the snakemake workflow work until the real function is ready
    It changes arbitrarily chosen values depending on slag content, not based on physics or anything.

    Parameters
    ----------
    slag_ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    phi_paste: ??
        input from Atuls parameter identification, not sure currently

    Returns
    -------
    paste_youngs_modulus : float / pint stress unit, will be in 'GPa'
        approximated youngs modulus of paste
    paste_strength : float / pint stress unit, will be in 'MPa'
        approximated compressive strength of paste
    """
    # TODO: now the phi inputs are dummy, need to add a fucntin there maybe.

    paste_youngs_modulus_min = 30 * ureg('GPa')
    paste_youngs_modulus_max = 60 * ureg('GPa')
    #paste_youngs_modulus = paste_youngs_modulus_min + (paste_youngs_modulus_max - paste_youngs_modulus_min) * slag_ratio
    # E is "mostly" inversely proportional to the slag
    #paste_youngs_modulus = paste_youngs_modulus_min + (paste_youngs_modulus_max - paste_youngs_modulus_min) * slag_ratio
    paste_youngs_modulus = paste_youngs_modulus_max - (paste_youngs_modulus_max - paste_youngs_modulus_min) * slag_ratio

    paste_strength_min = 25 * ureg('MPa')
    paste_strength_max = 40 * ureg('MPa')
    #paste_strength = paste_strength_min + (paste_strength_max - paste_strength_min) * slag_ratio
    # E is "mostly" inversely proportional to the slag
    #paste_strength = paste_strength_min + (paste_strength_max - paste_strength_min) * slag_ratio
    paste_strength = paste_strength_max - (paste_strength_max - paste_strength_min) * slag_ratio

    # ATUL : temporary q(b|x)~N(mu,cov), b=(sigma_paste,E_paste), x = slag ratio
    th.manual_seed(seed=seed)
    mean = [float(paste_youngs_modulus.magnitude), float(paste_strength.magnitude)]
    dist = th.distributions.MultivariateNormal(loc=th.as_tensor(mean), covariance_matrix=th.tensor(phi_cov))
    paste_youngs_modulus, paste_strength = dist.sample()
    return paste_youngs_modulus.item() * ureg('GPa') , paste_strength.item() * ureg('MPa')


if __name__ == "__main__":
    # test while developing this
    E, fc = dummy_paste_strength_stiffness(0.8, phi_mean=[[1., 25], [0., 1.]], phi_cov=[[1., 0], [0., 1.]], seed=5)
    print(E, fc)
