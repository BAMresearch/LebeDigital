import pytest

from lebedigital.demonstrator_scripts import beam_design
from lebedigital.unit_registry import ureg


def test_beam_design_value():
    """
    Test function to test beam design module to make sure the result don't change
    This example is analogus to one described in this blog:
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """

    width, height = beam_design.section_dimension_rule_of_thumb(span=6.75 * ureg("m"))
    results = beam_design.check_beam_design(
        span=6750 * ureg("mm"),
        width=width,
        height=height,
        point_load=36e3 * ureg("N"),
        distributed_load=0 * ureg("N/mm"),
        compr_str_concrete=20 * ureg("N/mm^2"),
        yield_str_steel=500 * ureg("N/mm^2"),
        steel_dia_bu=12 * ureg("mm"),
        cover_min=2.5 * ureg("cm"),
    )

    assert results["n_steel_bars"] == 5 * ureg("")
    assert results["diameter"] == 10 * ureg("mm")
    assert results["crosssection"].magnitude == pytest.approx(392.699082)
    assert results["crosssection"].units == ureg("mm^2")


def test_beam_design_constraint():
    results = beam_design.check_beam_design(
        span=6750 * ureg("mm"),
        width=150 * ureg("mm"),
        height=455 * ureg("mm"),
        point_load=36e3 * ureg("N"),
        distributed_load=0 * ureg("N/mm"),
        compr_str_concrete=5 * ureg("N/mm^2"),
        yield_str_steel=500 * ureg("N/mm^2"),
        steel_dia_bu=12 * ureg("mm"),
        cover_min=5 * ureg("cm"),
    )

    assert results["constraint_min_fc"].magnitude > 0.0
    assert results["constraint_max_steel_area"].magnitude > 0.0


def test_beam_design_cover():
    """
    Test function to test beam design module to make sure the result don't change
    This example is analogus to one described in this blog:
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """

    width, height = beam_design.section_dimension_rule_of_thumb(span=6.75 * ureg("m"))
    results_1 = beam_design.check_beam_design(
        span=6750 * ureg("mm"),
        width=width,
        height=height,
        point_load=36e3 * ureg("N"),
        distributed_load=0 * ureg("N/mm"),
        compr_str_concrete=20 * ureg("N/mm^2"),
        yield_str_steel=500 * ureg("N/mm^2"),
        steel_dia_bu=12 * ureg("mm"),
        cover_min=5 * ureg("mm"),
    )

    results_2 = beam_design.check_beam_design(
        span=6750 * ureg("mm"),
        width=width,
        height=height,
        point_load=36e3 * ureg("N"),
        distributed_load=0 * ureg("N/mm"),
        compr_str_concrete=20 * ureg("N/mm^2"),
        yield_str_steel=500 * ureg("N/mm^2"),
        steel_dia_bu=12 * ureg("mm"),
        cover_min=8 * ureg("mm"),
    )

    assert results_1["n_steel_bars"] == 7 * ureg("")
    assert results_1["diameter"] == 8 * ureg("mm")
    assert results_1["crosssection"].magnitude == pytest.approx(351.85837720205683)
