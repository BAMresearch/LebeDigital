from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.approximate_tensile_strength import approximate_tensile_strength
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

def test_approximate_tensile_strength():
    compressive_strength = 30 * ureg('MPa')

    tensile_strength = approximate_tensile_strength(compressive_strength)

    assert_approx(tensile_strength, compressive_strength/10, rtol=0.001)