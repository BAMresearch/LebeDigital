from lebedigital.unit_registry import ureg
import pytest

#@ureg.wraps('MPa', ('MPa', 'kg/m^3'))
def computation_volume_content(input) :
    """
    This is the function to compute volume contents based on mass or volume ratios

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    input : dic / pint units required converted to matching units
        required:
            - 'density_cement' 
            - 'density_water'
            - 'wb_mass_ratio' : water to binder ratio water/(cement + slag (or cemII))
            - 'density_aggregates'
            - 'aggregates_volume_fraction'
            - 'density_sub'
            - 'sc_volume_fraction'
        optional (otherwise assumed to be zero)
             - 'plasticizer_content'
             - 'density_plasticizer'
        
    Returns
    -------
    output : dic / pint unit given
             - 'water_vol_fraction'
             - 'sub_vol_fraction'
             - 'cem_vol_fraction'
             - 'cem_vol_content'
             - 'sub_vol_content'
             - 'water_vol_content'
             - 'aggregates_vol_content'
             - 'density_paste'
    """

    # initialize output dictionary
    output = {}

    # optional paramter
    if 'plasticizer_content' not in input.keys():
        input['plasticizer_content'] = 0 * ureg('kg/m^3')
        input['density_plasticizer'] = 42 * ureg('kg/m^3')  # dummy value

    # converting to correct pint units / automatic check for pint input / check for required input values
    input['density_sub'].ito('kg/m^3')
    input['density_cement'].ito('kg/m^3')
    input['density_water'].ito('kg/m^3')
    input['density_plasticizer'].ito('kg/m^3')
    input['density_aggregates'].ito('kg/m^3')

    # compute binder density
    density_binder = input['sc_volume_fraction']*input['density_sub']+(1-input['sc_volume_fraction'])*input['density_cement']

    # compute volume ratio of water to binder (cement and slag)
    water_vol_fraction_to_binder = input['wb_mass_ratio']*density_binder/\
                         (input['density_water']+input['wb_mass_ratio']*density_binder)

    # volume fractions
    ## plasticizer
    output['plasticizer_vol_fraction'] = input['plasticizer_content']/input['density_plasticizer']
    ## water
    output['water_vol_fraction'] = (1 - input['aggregates_volume_fraction'])*water_vol_fraction_to_binder\
                                   - output['plasticizer_vol_fraction']
    ## binder
    binder_vol_fraction = (1 - input['aggregates_volume_fraction'])*(1 - water_vol_fraction_to_binder)

    output['sub_vol_fraction'] = binder_vol_fraction*input['sc_volume_fraction']
    output['cem_vol_fraction'] = binder_vol_fraction*(1-input['sc_volume_fraction'])

    # sanity check
    assert 1 == pytest.approx((output['plasticizer_vol_fraction'] +
                              output['water_vol_fraction'] +
                              output['sub_vol_fraction'] +
                              output['cem_vol_fraction'] +
                              input['aggregates_volume_fraction']).magnitude)

    # computation of volume contents
    output['cem_vol_content'] = output['cem_vol_fraction'] * input['density_cement']
    output['sub_vol_content'] = output['sub_vol_fraction'] * input['density_sub']
    output['water_vol_content'] = output['water_vol_fraction'] * input['density_water']
    output['aggregates_vol_content'] = input['aggregates_volume_fraction'] * input['density_aggregates']

    # paste density
    output['density_paste'] = (output['cem_vol_content'] + output['sub_vol_content'] +
                              output['water_vol_content'] + input['plasticizer_content']) /\
                              (1-input['aggregates_volume_fraction'])

    return output



def computation_ratios(input):

    # set correct units
    input['density_cement'].ito('kg/m^3')
    input['density_slag'].ito('kg/m^3')
    input['density_water'].ito('kg/m^3')
    input['density_plasticizer'].ito('kg/m^3')
    input['density_aggregates'].ito('kg/m^3')

    input['plasticizer_content'].ito('kg/m^3')
    input['cem_vol_content'].ito('kg/m^3')
    input['slag_vol_content'].ito('kg/m^3')
    input['water_vol_content'].ito('kg/m^3')
    input['aggregates_vol_content'].ito('kg/m^3')

    output = {}

    # compute water to binder mass ratio
    output['wb_mass_ratio'] = (input['water_vol_content']+input['plasticizer_content'])/\
                              (input['cem_vol_content']+input['slag_vol_content'])

    # compute slag to cement volume fraction
    output['sc_volume_fraction'] = input['slag_vol_content']/(input['slag_vol_content'] + input['cem_vol_content'])




    # output['aggregates_volume_fraction']
    # output['sc_volume_fraction']


    # compute ratios


    return output




















if __name__ == "__main__":
    # test while developing this
    print('test')
    input = {}
    # densities
    input['density_cement'] = 1000 * ureg('kg/m^3')
    input['density_sub'] = 1000 * ureg('kg/m^3')
    input['density_water'] = 1000 * ureg('kg/m^3')
    input['density_plasticizer'] = 1000 * ureg('kg/m^3')
    input['density_aggregates'] = 1000 * ureg('kg/m^3')

    input['wb_mass_ratio'] = 0.4
    input['sc_volume_fraction'] = 0.4
    input['aggregates_volume_fraction'] = 0.65
    input['plasticizer_content'] = 10 * ureg('kg/m^3') # optional


    output = computation_volume_content(input)
    for key in output:
        print(key,':', output[key])
    #print(output)