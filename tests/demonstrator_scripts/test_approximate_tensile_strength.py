from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.approximate_tensile_strength import approximate_tensile_strength
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest

def test_approximate_tensile_strength():

    compressive_strength = 30 * ureg('MPa')
    tensile_strength = approximate_tensile_strength(compressive_strength)
    assert tensile_strength.magnitude == pytest.approx(2.896468153816889)

    compressive_strength = 60 * ureg('MPa')
    tensile_strength = approximate_tensile_strength(compressive_strength)
    assert tensile_strength.magnitude == pytest.approx(4.354742315434558)
