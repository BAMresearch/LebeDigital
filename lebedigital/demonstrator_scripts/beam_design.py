import math
from lebedigital.unit_registry import ureg
import numpy as np
import pint


@ureg.wraps(('mm', 'mm'), 'mm')
def section_dimension_rule_of_thumb(span: pint.Quantity) -> tuple[pint.Quantity, pint.Quantity]:
    """
    This is as per eurocode guideline choice of section dimension based on span
    Only used as a useful input for the test.

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    span :
        beam span in mm

    Returns
    -------
    tuple :
        tuple of width and height of the beam in mm.
    """
    # span/depth ratio for simply supported beam is 15
    height = 50 * round((span/15)/50)
    # width of beam
    width = 0.5 * height

    return width, height


@ureg.check('[length]', '[force]', '[force]/[length]')
def max_bending_moment_and_shear_force(span: pint.Quantity,
                                       point_load: pint.Quantity,
                                       distributed_load: pint.Quantity) -> tuple[pint.Quantity, pint.Quantity]:
    """
    function to compute max bending moment and shear force for simply supported beam with point load or distributed load

    The input is checked to be using the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    span :
        beam length
    point_load :
        point load
    distributed_load :
        distributed load

    Returns
    -------
    tuple :  pint force and moment units, specific units depending on input
        tuple of  maximum moment and maximum shear force
    """

    max_moment_dist_load = distributed_load * span**2 / 8
    max_shear_force_dist_load = distributed_load*span / 2
    max_moment_point_load = point_load*span / 4
    max_shear_force_point_load = point_load / 2

    return (max_moment_point_load + max_moment_dist_load, 
            max_shear_force_dist_load + max_shear_force_point_load)


@ureg.wraps(('mm^2',''), ('mm', 'mm', 'N*mm', 'N/mm^2', 'N/mm^2', 'mm', 'mm', 'mm'))
def beam_required_steel(
        width: pint.Quantity,
        height: pint.Quantity,
        max_moment: pint.Quantity,
        fck: pint.Quantity,
        fyk: pint.Quantity,
        steel_dia: pint.Quantity,
        steel_dia_bu: pint.Quantity,
        cover: pint.Quantity,
) -> pint.Quantity:
    """
    Function to design singly reinforced beam with minimum shear reinforcement required.

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    width: float / pint length unit, will be converted to 'mm'
        beam width in mm.
    height: float / pint length unit, will be converted to 'mm'
        beam depth in mm.
    max_moment : float / pint moment unit, will be converted to 'N*mm'
        Maximum bending moment in N-mm.
    fck : float / pint stress unit, will be converted to 'N/mm^2'
        charateristic compressive strength of concrete in N/mm2.
    fyk : float / pint stress unit, will be converted to 'N/mm^2'
        Yield strength of steel in N/mm2.
    steel_dia : int
        Diameter of steel in mm for longitudinal reinforcement.
    steel_dia_bu : int
        Diameter of steel in mm, for Buegel-Bewehrung.
    cover : float / pint length unit, will be converted to 'mm'
        Depth of cover required as per exposure class in mm.

    Returns
    -------
    req_steel : with fixed pint units when appropriate (length in 'mm')
        Design of the reinforced beam section.
    fc_constraint : checks if fc is fallen below the minimum. negative values are bad
    """
    # effective section depth
    deff = height - cover - steel_dia_bu - steel_dia / 2
    # fcd=Design compressive strength
    a_cc = 0.85
    gamma_c = 1.5  # Concrete partial material safety factor
    fcd = a_cc * fck / gamma_c  # N/mm^2
    gamma_s = 1.15
    fywd = fyk / gamma_s  # N/mm^2
    # Bending measurement (here with stress block) (Biegebemessung (hier mit Spannungsblock))
    mued = max_moment / (width * deff ** 2 * fcd)

    # when mued >= 0.5, xi cannot be computed, the compressive strength is to low
    # this causes problems for the optimization scheme
    # therefore we set an effective mued, which is wrong, but we have a constraint to check for that
    if mued >= 0.5:
        mued_eff = 0.5
    else:
        mued_eff = mued

    fc_constraint = 0.5 - mued

    xi = 0.5 * (1 + math.sqrt(1 - 2 * mued_eff))
    req_steel = 1 / fywd * max_moment / (xi * deff)

    return req_steel, fc_constraint


@ureg.check('[length]','[length]', '[length]','[length]')
def get_max_reinforcement(acceptable_reinforcement_diameters : list, width, cover_min, steel_dia_bu):
    """
    computes the maximum reinforcement, that fits, based on geometry

    Parameters
    ----------
    acceptable_reinforcement_diameters: list with possible diameters
    width : width of the beam
    cover_min : minimum concrete cover
    steel_dia_bu: diameter of stirrups

    Returns
    -------
    max_area: maximum area of steel reinforcement
    """
    for diameter in reversed(acceptable_reinforcement_diameters):
        if cover_min < diameter:
            cover = diameter
        else:
            cover = cover_min

        n_steel = 2 * ureg('')
        # if two bars are too much, go to the lower diameter
        if not beam_check_spacing(diameter, n_steel, steel_dia_bu, width, cover):
            continue

        # for a level, that fits at least two bars, see what the maximum is
        while beam_check_spacing(diameter, n_steel, steel_dia_bu, width, cover):
            max_area = n_steel * np.pi * (diameter / 2) ** 2
            n_steel += 1 * ureg('')

        break
    return max_area






@ureg.check('[length]', '[length]', '[length]', '[force]', '[force]/[length]',
            '[stress]', '[stress]', '[length]', '[length]')
def check_beam_design(span: pint.Quantity,
                      width: pint.Quantity,
                      height: pint.Quantity,
                      point_load: pint.Quantity,
                      distributed_load: pint.Quantity,
                      compr_str_concrete: pint.Quantity,
                      yield_str_steel: pint.Quantity,
                      steel_dia_bu: pint.Quantity,
                      cover_min: pint.Quantity,
                      ) -> dict[str, pint.Quantity]:
    """
    Function to check specified design for area of steel

    The input is checked to be using the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    span : float / pint length unit
        length of the beam in mm.
    width: float / pint length unit
        beam width in mm.
    height: float / pint length unit
        beam depth in mm.
    point_load: float / pint force unit
        point load in center of beam
    distributed_load: float / pint force/length unit
        constant load along the beam
    compr_str_concrete : float / pint stress unit
        charateristic compressive strength of concrete in N/mm2.
    yield_str_steel : float / pint stress unit
        Yield strength of steel in N/mm2.
    steel_dia_bu : float / pint length unit
        Diameter of steel in mm for Buegel-Bewehrung
    cover_min : float / pint length unit
        Depth of cover required as per exposure class in mm.

    Returns
    -------
    float
        normalized difference of  specified area and required area given the diameter of steel and number of steel bars
        in the bottom of the section. It is negative if  design is not satisfied and positive if design is satisfied.
        Optimal will be close to zero.
    """
    max_moment, max_shear_force = max_bending_moment_and_shear_force(span,
                                                                     point_load,
                                                                     distributed_load)

    discrete_reinforcement = {'crosssection': np.nan, 'n_steel_bars': np.nan, 'diameter': np.nan}

    acceptable_reinforcement_diameters = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 28.0, 32.0, 40.0] * ureg('mm')

    # get max reinforcement
    max_reinforcement = get_max_reinforcement(acceptable_reinforcement_diameters, width, cover_min, steel_dia_bu)

    for diameter in acceptable_reinforcement_diameters:

        # set correct cover
        if cover_min < diameter:
            cover = diameter
        else:
            cover = cover_min

        required_area, fc_error = beam_required_steel(width,
                                            height,
                                            max_moment,
                                            compr_str_concrete,
                                            yield_str_steel,
                                            diameter,
                                            steel_dia_bu,
                                            cover)

        area = (np.pi * (diameter / 2) ** 2)  # mm^2
        nsteel = max(2.0 * ureg(''), np.rint(required_area / area))  # rounds up

        # set constraints
        discrete_reinforcement['constraint_min_fc'] = fc_error
        discrete_reinforcement['constraint_max_steel_area'] = (max_reinforcement - required_area)/max_reinforcement

        # combined constraint
        if discrete_reinforcement['constraint_min_fc'] < 0.0 or\
           discrete_reinforcement['constraint_max_steel_area'] < 0.0:
            sign = -1
        else:
            sign = 1

        discrete_reinforcement['constraint_beam_design'] = (sign * abs(discrete_reinforcement['constraint_min_fc']) *
                                                            abs(discrete_reinforcement['constraint_max_steel_area']))

        if beam_check_spacing(diameter, nsteel, steel_dia_bu, width, cover):
            # found the smallest diameter that has correct spacing
            discrete_reinforcement['crosssection'] = area * nsteel
            discrete_reinforcement['n_steel_bars'] = nsteel
            discrete_reinforcement['diameter'] = diameter
            break
        else:
            discrete_reinforcement['crosssection'] = required_area
            discrete_reinforcement['n_steel_bars'] = 2
            discrete_reinforcement['diameter'] = 2 * np.sqrt(required_area / 2 / np.pi)

    return discrete_reinforcement


@ureg.check('[length]', None, '[length]', '[length]', '[length]')
def beam_check_spacing(diameter_l, n_steel, diameter_bu, width, cover) -> bool:
    """

    Parameters
    ----------
    diameter_l : diameter of the steel reinforcement in longitudinal direction
    n_steel : number of steel bars
    diameter_bu : diameter of the stirrups
    width : width of the beam
    cover : concrete cover

    Returns
    -------
    bool : True when there is space for the given reinforcement, False, when not
    """
    assert n_steel >= 2 * ureg('')
    # effective width for reinforcements
    b_eff = width - 2 * cover - 2 * diameter_bu
    # compute spacing
    s = (b_eff - n_steel * diameter_l)/(n_steel - 1)

    # set minimum spacing
    # currently ignoring aggregate size, diameter_largest_rock + 5mm is another constraint
    s_min = 20 * ureg('mm')
    if diameter_l > 20 * ureg('mm'):
        s_min = diameter_l

    if s_min > s:
        return False
    else:
        return True
