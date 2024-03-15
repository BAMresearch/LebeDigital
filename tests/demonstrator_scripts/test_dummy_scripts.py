import pytest

from lebedigital.demonstrator_scripts.dummy_hydration_parameters import \
    dummy_hydration_parameters
from lebedigital.demonstrator_scripts.dummy_paste_strength_stiffness import \
    dummy_paste_strength_stiffness


def test_dummy_scripts():
    # just to fix the coverage
    B1, B2, eta, E_act, Q_pot, T_ref = dummy_hydration_parameters(0.5, 10)
    assert B1.magnitude == pytest.approx(0.0002208)
    assert B2.magnitude == pytest.approx(0.0024229)
    assert eta.magnitude == pytest.approx(5.554)
    assert E_act.magnitude == pytest.approx(5653 * 8.3145)
    assert Q_pot.magnitude == pytest.approx(200000)
    assert T_ref.magnitude == pytest.approx(25)

    E, fc = dummy_paste_strength_stiffness(0, 10)
    assert E.magnitude == pytest.approx(60)
    assert fc.magnitude == pytest.approx(40)
