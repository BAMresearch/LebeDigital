
from lebedigital.BeamDesign import *

def test_beam_design():
    """
    Test function to test beam design module 
    """
    b,d = section_dimension_rule_of_thumb(span=6500) #mm
    MaxMoment,MaxShearForce =  max_bending_moment_and_shear_force(span=6500,load=30e3,load_type="point_load")
    out=  beam_section_design(span =6500,
                              b = b,
                              d = d,
                              MaxMoment=MaxMoment,
                              MaxShearForce=MaxShearForce,
                              fck=20,
                              fyk=500,
                              steelDia=12,
                              cover=25)
    for key,val in out.items():
        print(f"{key} --> {val}")

if __name__ =="__main__":
	test_beam_design()