import numpy as np
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
from lebedigital.demonstrator_scripts.youngs_modulus_approximation import youngs_modulus_approximation
from lebedigital.demonstrator_scripts.youngs_modulus_approximation import ureg

def test_youngs_modulus_approximation() :
    """
    This is the test for the youngs_modulus_aprroximation function
    """

    # regular test
    fc = 30 * ureg('MPa')
    density = 2400 * ureg('kg/m^3') # kg/m
    E = youngs_modulus_approximation(fc,density)

    assert_approx(E, 25439.083860411065* ureg('MPa'), rtol=0.001)

    # test with numpy array
    fc = np.array([30,90]) * ureg('MPa')
    density = ([2400,2400]) * ureg('kg/m^3')  # kg/m
    E = youngs_modulus_approximation(fc, density)

    assert_approx(E, np.array([25439.083860411065, 38750.98044651661]) * ureg('MPa'), rtol=0.001)

    # test with changed units
    fc = 30000000 * ureg('Pa')
    density = 2400 * ureg('kg/m^3') # kg/m
    E = youngs_modulus_approximation(fc,density)

    assert_approx(E, 25439.083860411065* ureg('MPa'), rtol=0.001)