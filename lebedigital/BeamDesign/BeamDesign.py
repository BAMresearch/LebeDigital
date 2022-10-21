# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 11:13:33 2022

@author: zv1174
"""
import math

def section_dimension_rule_of_thumb(span:float)->tuple:
    """
    This is as per eurocode guideline choice of section dimension based on span 

    Parameters
    ----------
    span : float
        beam span in mm.

    Returns
    -------
    tuple
        tuple of width and depth of the beam in mm.

    """
    #span/depth ratio for simply supported beam is 15 Therefore,
    d = 50*round((span/15)/50)
    #width of beam
    b = 0.5*d
    return b,d 

def max_bending_moment_and_shear_force(span:float,load:float,load_type:str) ->tuple:
    """
    function to compute max bending moment and shear force for simply supported beam with point load or distributed load

    Parameters
    ----------
    span : float
        beam span in mm.
    load : float
        beam load in N.
    load_type : str
        whether point_load or distributed_load to describe type of load.

    Returns
    -------
    tuple
        tuple of  maximum moment and maximum shear force.
    """
    w,l  = load,span
    if load_type=="distributed_load":
        max_moment=  w*l**2/8 
        max_shear_force = w*l/2
        return max_moment, max_shear_force
    if load_type=="point_load":
        max_moment=  w*l/4
        max_shear_force = w/2
        return max_moment, max_shear_force

def beam_section_design(
                       span:float,
                       b:float,
                       d:float,
                       MaxMoment:float,
                       MaxShearForce:float,
                       fck:float,
                       fyk:float,
                       steelDia:float,
                       cover:float,
                        ) -> dict:
    """
    Function to design singly reinforced beam with minimum shear reinforcement required.

    Parameters
    ----------
    span : float
        Span of the beam in mm.
    b: float
        beam width in mm.
    d: float
        beam depth in mm. 
    MaxMoment : float
        Maximum bending moment in N-mm.
    MaxShearForce : float
        Maximum shear force in N.
    fck : float
        charateristic compressive strength of concrete in N/mm2.
    fyk : float
        Yield strength of steel in N/mm2.
    steelDia : float
        Diameter of steel in mm.
    cover : float
        Depth of cover required as per exposure class in mm.

    Returns
    -------
    dict
        Design of the reinforced beam section.

    """
    #effective section depth
    deff = d - cover - steelDia - steelDia/2
    #fcd=Design compressive strength
    a_cc=0.85
    gamma_c=1.5 #Concrete partial material safety factor
    fcd=a_cc*fck/gamma_c #N/mm^2
    gamma_s=1.0
    fywd=fyk/gamma_s#N/mm^2
    #Bending measurement (here with stress block) (Biegebemessung (hier mit Spannungsblock))
    muEd=MaxMoment/(b*deff**2*fcd)/1000 
    xi=0.5*(1+math.sqrt(1-2*muEd))
    As1=1/fywd*MaxMoment/(xi*deff) #[-] 
    A=(math.pi*steelDia**2/4) #mm^2
    nsteel= math.ceil(As1/A) #rounds up
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
    bw=b #Rectangular cross section [mm]
    V_Rdmax= alpha_cw*bw*z*v1*fcd/(cot+1/cot) #[N]
    assert V_Rdmax > MaxShearForce, "Compression strut not stable"
    #Tension strut----
    #VEdred= Shear force reduction -> only possible is "direct support" 
    VEdred=MaxShearForce #here assumed: indirect support  [N]
    #a_sw =Steal stection per m 
    a_sw=VEdred/(z*fywd*cot)*1000 #[mm^2/m]
    #Calculated minimum shear force reinforcement----
    alpha=math.radians(90)# Angle between shear force reinforcement and the component axis perpendicular to the shear force: α = 90° (vertical stirrups)
    #fctm mean tensile strength 
    if fck<=50: fctm=0.3*fck**(2/3)
    else: fctm=2.12*math.ln(1+((fck+8)/10)) 
    a_swmin=0.16*fctm/fyk*(bw/1000)*math.sin(alpha)*1e6 #[mm2/m]
    #Statically required stirrup spacing----
    #print(f"stirrup and bending reinforcment have same diameter = {steelDia} [mm]")
    a_used= max(a_sw,a_swmin) #[mm^2/m]
    Aw=2*math.pi*12**2/4 #[mm²]
    slmax=Aw/a_used*1000 #[mm]
    #Minimum structural reinforcement----
    if fck<=50:
        if MaxShearForce<= 0.3*V_Rdmax:
            S= min(0.7*d, 300, slmax)
        elif MaxShearForce<= 0.6*V_Rdmax:
            S= min(0.5*d, 300, slmax)   
        else:
            S= min(0.25*d, 200, slmax)  
    else:
        if MaxShearForce<= 0.3*V_Rdmax:
            S= min(0.7*d, 200, slmax)
        elif MaxShearForce<= 0.6*V_Rdmax:
            S= min(0.7*d, 200, slmax)   
        else:
            S= min(0.25*d, 200, slmax)  
    #output final design
    out = {}
    out["Depth[mm]"]=d
    out["cover[mm]"] = cover
    out["width[mm]"] = b
    out["bottom_steel_dia[mm]"] = steelDia
    out["bottom_steel_numbers"] = nsteel
    out["top_steel_dia[mm]"] = steelDia
    out["top_steel_numbers"] = 2
    out["shear_reinforcement_dia[mm]"] = 12
    out["shear_reinforcement_spacing[mm]"] = S
    return out
       
