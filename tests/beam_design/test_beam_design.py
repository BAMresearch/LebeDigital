import sys,os
sys.path.append(os.path.join("..",".."))
from lebedigital.BeamDesign import BeamDesign as bd

def test_beam_design():
    """
    Test function to test beam design module 
    This example is analogus to one described in this blog: 
    https://www.structuralguide.com/worked-example-singly-reinforced-beam-design-using-ec2/
    """
    b,d = bd.section_dimension_rule_of_thumb(span=6750) #mm
    MaxMoment,MaxShearForce =  bd.max_bending_moment_and_shear_force(span=6750,load=36e3,load_type="point_load")
    out=  bd.beam_section_design(span =6750,
                              b = b,
                              d = d,
                              MaxMoment=MaxMoment,
                              MaxShearForce=MaxShearForce,
                              fck=20,
                              fyk=500,
                              steelDia=16,
                              cover=25)
    for key,val in out.items():
        print(f"{key} --> {val}")
    assert(out["bottom_steel_numbers"]==2)
    assert(out["shear_reinforcement_spacing[mm]"]==300)

if __name__ =="__main__":
	test_beam_design()
