from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.approximate_max_degree_of_hydration import approximate_max_degree_of_hydration
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest

def test_approximate_max_degree_of_hydration():

    water_cement_ratio = 0.0 * ureg('')
    max_degree_of_hydration = approximate_max_degree_of_hydration(water_cement_ratio)
    assert max_degree_of_hydration == pytest.approx(0.0)

    water_cement_ratio = 0.2 * ureg('')
    max_degree_of_hydration = approximate_max_degree_of_hydration(water_cement_ratio)
    assert max_degree_of_hydration == pytest.approx(0.523350254)

    water_cement_ratio = 0.4 * ureg('')
    max_degree_of_hydration = approximate_max_degree_of_hydration(water_cement_ratio)
    assert max_degree_of_hydration == pytest.approx(0.6942760942760942)

    water_cement_ratio = 1.0 * ureg('')
    max_degree_of_hydration = approximate_max_degree_of_hydration(water_cement_ratio)
    assert max_degree_of_hydration == pytest.approx(0.863484087)
