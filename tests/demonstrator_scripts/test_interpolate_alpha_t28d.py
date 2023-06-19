import pytest

from lebedigital.demonstrator_scripts.interpolate_alpha_t28d import interpolate_alpha_t28d


def test_interpolate_alpha_t28d():
    alpha_mix1 = 0.5
    alpha_mix2 = 1.0

    slag_ratio = 0.0
    assert interpolate_alpha_t28d(alpha_mix1, alpha_mix2, slag_ratio) == pytest.approx(alpha_mix1)
    slag_ratio = 0.5
    assert interpolate_alpha_t28d(alpha_mix1, alpha_mix2, slag_ratio) == pytest.approx((alpha_mix1 + alpha_mix2) / 2)
    slag_ratio = 1.0
    assert interpolate_alpha_t28d(alpha_mix1, alpha_mix2, slag_ratio) == pytest.approx(alpha_mix2)
