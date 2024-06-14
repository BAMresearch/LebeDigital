import pytest

from lebedigital.demonstrator_scripts.dummy_hydration_parameters import \
    dummy_hydration_parameters
from lebedigital.demonstrator_scripts.dummy_paste_strength_stiffness import \
    dummy_paste_strength_stiffness


def test_dummy_scripts():
    # just to fix the coverage
    phi_mean = [[0.01 * 2.916E-4, 2.916E-4], [0.01 * 0.0024229, 0.0024229], [0.01 * 5.554, 5.554],
                [0.01 * 5653 * 8.3145, 5653 * 8.3145], [-100000, 300000], [0.01 * 25, 25]]
    #phi_cov = np.diag(0.0025 * np.array(phi_mean)[:, 1]).tolist()
    phi_cov = [
            [
                7.29e-10,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0
            ],
            [
                0.0,
                6.05725e-08,
                0.0,
                0.0,
                0.0,
                0.0
            ],
            [
                0.0,
                0.0,
                0.013885000000000002,
                0.0,
                0.0,
                0.0
            ],
            [
                0.0,
                0.0,
                0.0,
                117.50467125000002,
                0.0,
                0.0
            ],
            [
                0.0,
                0.0,
                0.0,
                0.0,
                750.0,
                0.0
            ],
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0625
            ]
        ]
    seed = 10
    B1, B2, eta, E_act, Q_pot, T_ref = dummy_hydration_parameters(slag_ratio=0.5, phi_mean=phi_mean, phi_cov=phi_cov, seed=seed)
    assert B1.magnitude == pytest.approx(0.0002768,rel=1e-2)
    assert B2.magnitude == pytest.approx(0.002185,rel=1e-2)
    assert eta.magnitude == pytest.approx(5.5461,rel=1e-2)
    assert E_act.magnitude == pytest.approx(47223.5698,rel=1e-2)
    assert Q_pot.magnitude == pytest.approx(250025.191,rel=1e-2)
    assert T_ref.magnitude == pytest.approx(25.035,rel=1e-2)

    E, fc = dummy_paste_strength_stiffness(slag_ratio = 0, phi_mean=[[1., 25], [0., 1.]], phi_cov=[[1., 0], [0., 1.]], seed=5)
    assert E.magnitude == pytest.approx(59.51,rel=1e-2)
    assert fc.magnitude == pytest.approx(39.39,rel=1e-2)
