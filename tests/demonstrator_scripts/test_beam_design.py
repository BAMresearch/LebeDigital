from lebedigital.demonstrator_scripts import beam_design
import pytest
from lebedigital.unit_registry import ureg


@pytest.mark.parametrize("compr_str_concrete", [20, 60])
@pytest.mark.parametrize("point_load", [36e3, 60e3])
def test_beam_design_function(compr_str_concrete,point_load):
    """
    Test function to test beam design module to check if it runs
    """

    span = 6.75 * ureg('m')
    width, height = beam_design.section_dimension_rule_of_thumb(span=span)
    out = beam_design.check_design(span=span,
                                   width=width,
                                   height=height,
                                   point_load = point_load*ureg('N'),
                                   distributed_load= 0*ureg('N/mm'),
                                   compr_str_concrete=compr_str_concrete*ureg('N/mm^2'),
                                   yield_str_steel=500*ureg('N/mm^2'),
                                   steel_dia=12*ureg('mm'),
                                   n_bottom=8,
                                   cover=2.5*ureg('cm'))

    assert(out > 0)


def test_beam_design_value():
    """
    Test function to test beam design module to make sure the result don't change
    This example is analogus to one described in this blog: 
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """

    width, height = beam_design.section_dimension_rule_of_thumb(span=6.75*ureg('m'))
    out = beam_design.check_design(span=6750*ureg('mm'),
                                   width=width,
                                   height=height,
                                   point_load = 36e3*ureg('N'),
                                   distributed_load= 0*ureg('N/mm'),
                                   compr_str_concrete=20*ureg('N/mm^2'),
                                   yield_str_steel=500*ureg('N/mm^2'),
                                   steel_dia=12*ureg('mm'),
                                   n_bottom=4,
                                   cover=2.5*ureg('cm'))

    assert(out == pytest.approx(0.214974352))
