import pytest

from lebedigital.demonstrator_scripts.computation_loads_with_safety import computation_loads_with_safety
from lebedigital.unit_registry import ureg


def test_computation_loads_with_safety():
    safety_factor_permanent = 1.35 * ureg("")
    safety_factor_variable = 1.5 * ureg("")
    distributed_load_permanent = 15 * ureg("kN/m")
    distributed_load_variable = 20 * ureg("kN/m")
    point_load_permanent = 1 * ureg("N")
    point_load_variable = 1 * ureg("N")

    distributed_load, point_load = computation_loads_with_safety(
        safety_factor_permanent,
        safety_factor_variable,
        distributed_load_permanent,
        distributed_load_variable,
        point_load_permanent,
        point_load_variable,
    )

    assert pytest.approx(distributed_load.magnitude) == 50.25
    assert pytest.approx(point_load.magnitude) == 2.85
