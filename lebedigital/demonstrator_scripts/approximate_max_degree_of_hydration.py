from lebedigital.unit_registry import ureg
import numpy as np

def approximate_max_degree_of_hydration(wc):
    """
    function to approximate maximum degree of hydration
    the equation is published in Mills (1966): Factors influencing cessation of hydration in water cured cement pastes

    Parameters
    ----------
    wc : float
        water to cement ratio

    Returns
    -------
    alpha_max : float
        maximum degree of hydration
    """

    alpha_max = 1.031 * wc / (0.194 + wc)

    return alpha_max