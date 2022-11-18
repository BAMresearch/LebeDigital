import sys,os
sys.path.append(os.path.join("..",".."))
from lebedigital.BeamDesign import BeamDesign as bd
import pytest
def test_beam_design():
    """
    Test function to test beam design module 
    This example is analogus to one described in this blog: 
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """
    b,d = bd.section_dimension_rule_of_thumb(span=6750) #mm
    out=  bd.check_design(span =6750,
                              b = b,
                              d = d,
                              point_load = 36e3,
                              distributed_load= 0,
                              fck=20,
                              fyk=500,
                              steelDia=12,
                              n_bottom=3,
                              cover=25)
    assert(out == pytest.approx(40.7447380827))
if __name__ =="__main__":
	test_beam_design()
