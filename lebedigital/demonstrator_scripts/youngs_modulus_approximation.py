import numpy as np
from lebedigital.unit_registry import ureg

@ureg.wraps('MPa', ('MPa', 'kg/m^3'))
def youngs_modulus_approximation(fc,density) :
    """
    This is the function to approximate Young's modulus from compressive strength, according to ACI 363

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    This source was chosen based on the paper:
    Sun, Fanourakis: The validation of elastic modulus models: Code models and their modified versions, 2021
    DOI: 10.1002/suco.202100312

    Parameters
    ----------
    fc : float / pint stress unit, will be converted to 'MPa'
        concrete compressive strength after 28 days
    density: float / pint density unit, will be converted to 'kg/m^3'
        concrete density

    Returns
    -------
    youngs_modulus : float / pint stress unit, will be in 'MPa'
        approximated youngs modulus
    """

    youngs_modulus = 3320*np.sqrt(fc)+6895*(density/2320)**1.5

    return youngs_modulus


if __name__ == "__main__":
    # test while developing this
    fc = np.array([30,90]) * ureg('MPa')
    density = np.array([2400,2400]) * ureg('kg/m^3') # kg/m

    E = youngs_modulus_approximation(fc,density)
    print(E)
