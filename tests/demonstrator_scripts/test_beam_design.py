from lebedigital.demonstrator_scripts import beam_design
import pytest


def test_beam_design():
    """
    Test function to test beam design module 
    This example is analogus to one described in this blog: 
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """

    width, height = beam_design.section_dimension_rule_of_thumb(span=6750) # mm
    out = beam_design.check_design(span=6750,
                                   width=width,
                                   height=height,
                                   point_load = 36e3,
                                   distributed_load= 0,
                                   compr_str_concrete=20,
                                   yield_str_steel=500,
                                   steel_dia=12,
                                   n_bottom=3,
                                   cover=25)

    assert(out == pytest.approx(0.13647667348210732))


if __name__ == "__main__":
    test_beam_design()
