from lebedigital.demonstrator_scripts import beam_design
import pytest
from lebedigital.unit_registry import ureg


def test_beam_design():
    """
    Test function to test beam design module 
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
                                   n_bottom=3,
                                   cover=2.5*ureg('cm'))

    assert(out == pytest.approx(0.0479153789))
