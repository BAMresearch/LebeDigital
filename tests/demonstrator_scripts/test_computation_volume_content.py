import numpy as np
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
from lebedigital.demonstrator_scripts.computation_volume_content import computation_volume_content
from lebedigital.demonstrator_scripts.computation_volume_content import computation_ratios
from lebedigital.unit_registry import ureg

def test_computation_volume_content() :
    """
    This is the test for the computation_volume_content function.
    It tests a number of sanity checks as well as different input options
    """

    #########################################
    ### OPTION 1 : with slag to cem ratio ###
    #########################################
    # TEST No 1, check intended use in general
    input1 = {}
    # densities
    input1['density_cement'] = 1.44 * ureg('g/cm^3')
    input1['density_slag'] = 840 * ureg('kg/m^3')
    input1['density_water'] = 977 * ureg('kg/m^3')
    input1['density_plasticizer'] = 0.98 * ureg('g/cm^3')
    input1['density_aggregates'] = 1.5 * ureg('t/m^3')

    input1['wb_mass_ratio'] = 0.4
    input1['sc_volume_fraction'] = 0.4
    input1['aggregates_volume_fraction'] = 0.65
    input1['plasticizer_content'] = 10 * ureg('kg/m^3')

    output1 = computation_volume_content(input1)

    # check all values
    assert_approx(output1['water_vol_fraction'], 0.11960918301641323, rtol=0.001)
    assert_approx(output1['slag_vol_fraction'], 0.08807469414037347, rtol=0.001)
    assert_approx(output1['cem_vol_fraction'], 0.1321120412105602, rtol=0.001)
    assert_approx(output1['cem_vol_content'], 190.24133934320665 * ureg('kg/m^3'), rtol=0.001)
    assert_approx(output1['slag_vol_content'], 73.98274307791371 * ureg('kg/m^3'), rtol=0.001)
    assert_approx(output1['water_vol_content'], 116.85817180703572 * ureg('kg/m^3'), rtol=0.001)
    assert_approx(output1['aggregates_vol_content'], 975.0 * ureg('kg/m^3'), rtol=0.001)
    assert_approx(output1['density_paste'], 1117.3778692233032 * ureg('kg/m^3'), rtol=0.001)

    # TEST No 2, sanity check density paste
    input2 = {}
    # densities
    input2['density_cement'] = 1000 * ureg('kg/m^3')
    input2['density_slag'] = 1000 * ureg('kg/m^3')
    input2['density_water'] = 1000 * ureg('kg/m^3')
    input2['density_plasticizer'] = 1000 * ureg('kg/m^3')
    input2['density_aggregates'] = 1000 * ureg('kg/m^3')

    input2['wb_mass_ratio'] = 0.4
    input2['sc_volume_fraction'] = 0.4
    input2['aggregates_volume_fraction'] = 0.65
    input2['plasticizer_content'] = 10 * ureg('kg/m^3')

    output2 = computation_volume_content(input2)

    # check paste density
    assert_approx(output2['density_paste'], 1000 * ureg('kg/m^3'), rtol=0.001)


    # TEST No 3, sanity check w/b ratio of 1, compare weight with equal density
    input3 = {}
    # densities
    input3['density_cement'] = 1000 * ureg('kg/m^3')
    input3['density_slag'] = 840 * ureg('kg/m^3')
    input3['density_water'] = 1000 * ureg('kg/m^3')
    input3['density_plasticizer'] = 0.98 * ureg('g/cm^3')
    input3['density_aggregates'] = 1.5 * ureg('t/m^3')

    input3['wb_mass_ratio'] = 1
    input3['sc_volume_fraction'] = 0
    input3['aggregates_volume_fraction'] = 0

    output3 = computation_volume_content(input3)

    # check all values
    assert_approx(output3['water_vol_fraction'], output3['cem_vol_fraction'], rtol=0.001)

    #########################################
    ### OPTION 2 : with cemI to cemII ratio #
    #########################################
    # TEST No 4, check intended use in general
    input4 = {}
    # densities
    input4['density_cement'] = 1.0 * ureg('g/cm^3')
    input4['density_cemII'] = 2.0 * ureg('g/cm^3')
    input4['density_water'] = 977 * ureg('kg/m^3')
    input4['density_plasticizer'] = 0.98 * ureg('g/cm^3')
    input4['density_aggregates'] = 1.5 * ureg('t/m^3')

    input4['wb_mass_ratio'] = 0.4
    input4['cem_volume_ratio'] = 0.5
    input4['aggregates_volume_fraction'] = 0.65
    input4['plasticizer_content'] = 10 * ureg('kg/m^3')

    output4 = computation_volume_content(input4)

    # check cemI and cemII values
    assert_approx(output4['cem_vol_content']*2,output4['cemII_vol_content'], rtol=0.001)





def test_computation_ratio():
    # TEST No 1, check intended use in general
    input1 = {}
    # densities
    input1['density_cement'] = 1500 * ureg('kg/m^3')
    input1['density_slag'] = 500 * ureg('kg/m^3')
    input1['density_water'] = 1000 * ureg('kg/m^3')
    input1['density_plasticizer'] = 0.98 * ureg('g/cm^3')
    input1['density_aggregates'] = 1.5 * ureg('t/m^3')

    input1['wb_mass_ratio'] = 1
    input1['sc_volume_fraction'] = 0.5
    input1['aggregates_volume_fraction'] = 0.0
    input1['plasticizer_content'] = 0 * ureg('kg/m^3')

    output1 = computation_volume_content(input1)

    print('___________')
    print('___________')
    print('water_content : ', output1['water_vol_content'])
    print('cement_content : ', output1['cem_vol_content'])
    print('slag_content : ', output1['slag_vol_content'])
    print('___________')



    # TEST No 1, check intended use in general
    input2 = {}
    # densities
    input2['density_cement'] = input1['density_cement']
    input2['density_slag'] = input1['density_slag']
    input2['density_water'] =  input1['density_water']
    input2['density_plasticizer'] = input1['density_plasticizer']
    input2['density_aggregates'] = input1['density_aggregates']

    input2['plasticizer_content'] = input1['plasticizer_content']
    input2['cem_vol_content'] = output1['cem_vol_content']
    input2['slag_vol_content'] = output1['slag_vol_content']
    input2['water_vol_content'] = output1['water_vol_content']
    input2['aggregates_vol_content'] = output1['aggregates_vol_content']

    output2 = computation_ratios(input2)
    print('wb : ',input1['wb_mass_ratio'], output2['wb_mass_ratio'], )
    print('sc : ',input1['sc_volume_fraction'], output2['sc_volume_fraction'])
    print('___________')

    assert_approx(input1['wb_mass_ratio'], output1['wb_mass_ratio'], rtol=0.001)
    assert_approx(input1['sc_volume_fraction'], output1['sc_volume_fraction'], rtol=0.001)







# 
# 
# 
# 
# 
# def test_computation_ratio():
# 
#     # TEST No 1, check intended use in general
#     input1 = {}
#     # densities
#     input1['density_cement'] = 1000 * ureg('kg/m^3')
#     input1['density_slag'] = 1000 * ureg('kg/m^3')
#     input1['density_water'] = 1000 * ureg('kg/m^3')
#     input1['density_plasticizer'] = 0.98 * ureg('g/cm^3')
#     input1['density_aggregates'] = 1.5 * ureg('t/m^3')
# 
#     input1['plasticizer_content'] = 0 * ureg('kg/m^3')
#     input1['cem_vol_content'] = 100 * ureg('kg/m^3')
#     input1['slag_vol_content'] = 100 * ureg('kg/m^3')
#     input1['water_vol_content'] = 200 * ureg('kg/m^3')
#     input1['aggregates_vol_content'] = 975.0 * ureg('kg/m^3')
# 
#     output1 = computation_ratios(input1)
# 
#     print(output1)
# 
#     # TODO: FIX THIS TEST!!!!!!
# 
#     assert_approx(output1['wb_mass_ratio'], 0.4, rtol=0.001)
