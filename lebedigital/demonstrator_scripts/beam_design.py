import math
from lebedigital.unit_registry import ureg

@ureg.wraps(('mm','mm'), 'mm')
def section_dimension_rule_of_thumb(span: float) -> tuple:
    """
    This is as per eurocode guideline choice of section dimension based on span
    Only used as a useful input for the test.

    The input and output is wrapped by the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    span : float / pint length unit, will be converted to 'mm'
        beam span in mm

    Returns
    -------
    tuple : float / pint length unit in 'mm'
        tuple of width and height of the beam in mm.
    """
    #span/depth ratio for simply supported beam is 15 Therefore,
    height = 50 * round((span/15)/50)
    #width of beam
    width = 0.5 * height

    return width, height


@ureg.check('[length]','[force]','[force]/[length]')
def max_bending_moment_and_shear_force(span: float, point_load: float, distributed_load: float) -> tuple:
    """
    function to compute max bending moment and shear force for simply supported beam with point load or distributed load

    The input is checked to be using the python pint package: https://pint.readthedocs.io/
    This requires units to be attached to the input values.

    Parameters
    ----------
    span : float / pint length unit
        beam length in mm
    point_load : float / pint force unit
        point load in N
    distributed_load : str / pint force/length unit
        distributed load in N/mm

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


@ureg.wraps(None, ('mm', 'mm', 'N*mm', 'N', 'N/mm^2','N/mm^2','mm','mm'))
def beam_section_design(
                       width: float,
                       height: float,
                       max_moment: float,
                       max_shear_force: float,
                       fck: float,
                       fyk: float,
                       steel_dia: float,
                       cover: float,
                       ) -> dict:
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
    max_shear_force : float / pint force unit, will be converted to 'mm'
        Maximum shear force in N.
    fck : float / pint stress unit, will be converted to 'N/mm^2'
        charateristic compressive strength of concrete in N/mm2.
    fyk : float / pint stress unit, will be converted to 'N/mm^2'
        Yield strength of steel in N/mm2.
    steel_dia : int
        Diameter of steel in mm.
    cover : float / pint length unit, will be converted to 'mm'
        Depth of cover required as per exposure class in mm.

    Returns
    -------
    dict : with fixed pint units when appropriate (length in 'mm')
        Design of the reinforced beam section.
    """
    #effective section depth
    # TODO: first steel diameter is reinforcement in width direction. does it need to be the same? do we need to account for it at all?
    deff = height - cover - steel_dia - steel_dia / 2
    #fcd=Design compressive strength
    a_cc=0.85
    gamma_c=1.5 #Concrete partial material safety factor
    fcd=a_cc*fck/gamma_c #N/mm^2
    gamma_s=1.15
    fywd=fyk/gamma_s#N/mm^2
    #Bending measurement (here with stress block) (Biegebemessung (hier mit Spannungsblock))
    muEd= max_moment / (width * deff ** 2 * fcd)
    xi=0.5*(1+math.sqrt(1-2*muEd))
    As1= 1 / fywd * max_moment / (xi * deff) #[-]
    A=(math.pi * steel_dia ** 2 / 4) #mm^2
    nsteel= max(2, math.ceil(As1/A)) #rounds up
    #Compression strut angle pure bending (Druckstrebenwirkel)
    cot=1.2 #[-] cot(40)=1.2
    #Compression strut----
    #v1 = Reduction coefficient for concrete strength in shear cracks 
    if fck > 50: v1=1.0 
    else: v1=0.75
    #alpha_cw= Coefficient to account for the stress state in the compression chord (according to NA: 𝛼𝑐𝑤=1.0);
    alpha_cw=1.0
    # z=inner lever arm (innerer Hebelarm)
    c_vl=cover
    z=min(0.9*deff, max(deff-c_vl-30, deff-2*c_vl)) #in [mm]
    bw=width #Rectangular cross section [mm]
    V_Rdmax= alpha_cw*bw*z*v1*fcd/(cot+1/cot) #[N]
    assert V_Rdmax > max_shear_force, "Compression strut not stable"
    #Tension strut----
    #VEdred= Shear force reduction -> only possible is "direct support" 
    VEdred=max_shear_force #here assumed: indirect support  [N]
    #a_sw =Steal stection per m 
    a_sw=VEdred/(z*fywd*cot)*1000 #[mm^2/m]
    #Calculated minimum shear force reinforcement----
    alpha=math.radians(90)# Angle between shear force reinforcement and the component axis perpendicular to the shear force: α = 90° (vertical stirrups)
    #fctm mean tensile strength 
    if fck<=50: fctm=0.3*fck**(2/3)
    else: fctm=2.12*math.log(1+((fck+8)/10))
    a_swmin=0.16*fctm/fyk*(bw/1000)*math.sin(alpha)*1e6 #[mm2/m]
    #Statically required stirrup spacing----
    #print(f"stirrup and bending reinforcment have same diameter = {steelDia} [mm]")
    a_used= max(a_sw,a_swmin) #[mm^2/m]
    Aw=2*math.pi*12**2/4 #[mm²]
    slmax=Aw/a_used*1000 #[mm]
    #Minimum structural reinforcement----
    if fck<=50:
        if max_shear_force<= 0.3*V_Rdmax:
            S= min(0.7 * height, 300, slmax)
        elif max_shear_force<= 0.6*V_Rdmax:
            S= min(0.5 * height, 300, slmax)
        else:
            S= min(0.25 * height, 200, slmax)
    else:
        if max_shear_force<= 0.3*V_Rdmax:
            S= min(0.7 * height, 200, slmax)
        elif max_shear_force<= 0.6*V_Rdmax:
            S= min(0.7 * height, 200, slmax)
        else:
            S= min(0.25 * height, 200, slmax)
    #output final design
    out = {}
    out["top_steel_dia"] = steel_dia * ureg('mm')
    out["top_steel_numbers"] = 2
    out["shear_reinforcement_dia"] = 12 * ureg('mm')
    out["shear_reinforcement_spacing"] = S * ureg('mm')
    out["required_area_bottom_steel"] = As1 * ureg('mm^2')

    return out


@ureg.check('[length]','[length]','[length]','[force]','[force]/[length]','[stress]','[stress]','[length]',None,'[length]')
def check_design(span:float,
                 width:float,
                 height:float,
                 point_load:float,
                 distributed_load:float,
                 compr_str_concrete:float,
                 yield_str_steel:float,
                 steel_dia:float,
                 n_bottom:int,
                 cover:float,
                 ) -> float:
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
    steel_dia : float / pint length unit
        Diameter of steel in mm.
    n_bottom:int,
        Number of steel bars in the bottom of the section.
        On top for singly reinforced section we assume 2 steel bars of 12 mm which holds the cage
    cover : float / pint length unit
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

    design = beam_section_design(width,
                                 height,
                                 max_moment,
                                 max_shear_force,
                                 compr_str_concrete,
                                 yield_str_steel,
                                 steel_dia ,
                                 cover)

    specified_area = (math.pi * (steel_dia) ** 2 / 4) * n_bottom
    required_area = design["required_area_bottom_steel"]

    return (specified_area-required_area)/required_area
