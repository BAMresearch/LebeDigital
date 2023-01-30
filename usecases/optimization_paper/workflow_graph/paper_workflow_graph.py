import graphviz

import matplotlib.pyplot as plt
import os
#from pathlib import Path
import pathlib

# looks better but different to other text
plt.rc('mathtext', fontset='cm')  # matplotlib: force computer modern font set

LOCATION = pathlib.Path(__file__).parent.resolve()
IMAGE_FOLDER = LOCATION/'graph_tex_nodes'
pathlib.Path(IMAGE_FOLDER).mkdir(exist_ok=True)

FONTSIZE = 12


def tex2svg(file_name, formula, fontsize=14, dpi=300):
    """Render TeX formula to SVG and saves to file
       this uses the matplotlib tex functions:
       https://matplotlib.org/stable/tutorials/text/mathtext.html

    Args:
        file_name (str): file name
        formula (str): TeX formula.
        fontsize (int, optional): Font size.
        dpi (int, optional): DPI.
    Returns:
        full_file_name (str): path to image
    """

    full_file_name = IMAGE_FOLDER / (file_name + '.svg')

    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, formula, fontsize=fontsize)
    fig.savefig(full_file_name,dpi=dpi, transparent=True, format='svg',
                bbox_inches='tight', pad_inches=0.2)

    return full_file_name


def tex_node(dot_object,node_name,node_text, color=None, shape=None, fontsize=FONTSIZE):
    file_name = tex2svg(node_name,node_text,fontsize)
    dot_object.node(node_name, '',image=str(file_name), color=color,shape=shape)


def paper_workflow_graph(file_name = 'test_output', view=False):

    file_name = str(file_name) # convert pathlib object to useful string

    process = 'blue'
    input = 'green'
    problem = 'red'
    kpi = 'orange'

    dot = graphviz.Digraph(file_name, comment='Optimization Paper', format='pdf')
    dot.attr("node", fontsize=str(FONTSIZE))

    dot.node('epd products', 'environmental product decleration (epd) for cement, slag, aggregates in CO2/kg', color=input)
    dot.edge('epd products','co2 computation')
    dot.node('co2 computation', 'computation of global warming potential (GWP) per m³', color=process)
    dot.edge('co2 computation','gwp')
    dot.node('gwp', 'GWP per m³ [kg CO2,eq/m³]')
    dot.edge('gwp','co2 computation volume')
    dot.node('co2 computation volume', 'GWP per part', color=process)
    dot.edge('geometry','co2 computation volume')
    dot.edge('co2 computation volume','gwp total')
    dot.edge('cemI','co2 computation')
    dot.node('gwp total', 'GWP per part [kg CO2,eq]', color=kpi, shape='rectangle')
    dot.node('cemI', 'cement, slag, aggregate content \n[kg/m³]')
    dot.edge('volume computation', 'cemI')
    dot.edge('ratio_cemI_cemII','volume computation')
    dot.node('aggregate content', 'aggregate/(water+cement volume ratio) [-]', color=input, shape='rectangle')
    dot.node('ratio_cemI_cemII', 'slag to cement ratio', color=input, shape='rectangle')
    dot.edge('w/c','volume computation')
    dot.node('w/c', 'water cement mass ratio [-]', color=input)
    dot.node('plasticizer', 'plasticizer content [kg/m³]', color=input)
    dot.edge('plasticizer','volume computation')
    dot.edge('w/c','compute max doh')
    dot.node('compute max doh', 'function: max degree of hydration', color=process)
    dot.edge('compute max doh', 'max doh')
    dot.node('max doh', 'max degree of hydration [-]')
    dot.edge('max doh', 'fem model')
    dot.node('volume computation', 'computation of volume contents', color=process)
    dot.edge('aggregate content', 'concrete homogenization')
    dot.edge('ratio_cemI_cemII',  'interpolation_paste')
    dot.node('concrete strength 28d', 'concrete compressive strength \n28 days [N/mm2]')
    dot.edge('ft(fc)', 'concrete tensile strength 28d')
    dot.node('concrete tensile strength 28d', 'concrete tensile strength \n28 days [N/mm2]')
    dot.node('ft(fc)', 'function ft(fc)', color=process)
    dot.edge('concrete strength 28d', 'load bearing capacity')
    dot.edge('concrete strength 28d', 'ft(fc)')
    dot.edge('geometry', 'load bearing capacity')
    dot.node('load bearing capacity', 'beam design', color= process)
    dot.edge('load bearing capacity', 'kpi load bearing capacity')
    dot.node('kpi load bearing capacity', 'relative A diff:\n(A_prop-A_req)/A_req', color=kpi, shape='rectangle')
    dot.node('fem model', 'fem model', color=process)
    dot.edge('fem model','fem output')
    dot.node('fem output', 'max temperature over time\nmax yield over time',)
    dot.edge('fem output','fem postprocessing')
    dot.node('max temp', 'temperature constraint', color=input)
    dot.edge('max temp','fem postprocessing')
    dot.node('fem postprocessing', 'fem postprocessing', color=process)
    dot.edge('fem postprocessing','demolding time')
    dot.node('demolding time', 'time of demolding [h]', color=kpi, shape='rectangle')
    dot.edge('fem postprocessing','temperature check')
    dot.node('temperature check', 'temperature check [bool]', color=kpi, shape='rectangle')
    dot.node('concrete homogenization', 'concrete homogenization', color=process)
    dot.edge('concrete tensile strength 28d','fem model')
    dot.node('aggregate data', 'aggregate data: E, nu, C, kappa', color=input)
    dot.node('aggregate rho', 'aggregate density', color=input)
    dot.node('cem rho', 'densities: cement, slag, water', color=input)
    dot.edge('cem rho','volume computation')
    dot.edge('volume computation','paste rho')
    dot.node('paste rho', 'paste density')
    dot.edge('paste rho',  'concrete homogenization')
    dot.edge('aggregate data', 'concrete homogenization')
    dot.edge('aggregate rho',  'concrete homogenization')
    dot.edge('param paste strength 28d','concrete homogenization')
    dot.edge('aggregate rho', 'volume computation')
    dot.node('paste data', 'paste data: nu, C, kappa', color=input)
    dot.edge('paste data', 'concrete homogenization')

    dot.node('param paste strength 28d', 'paste strength and Youngs modulus\nafter 28 days [N/mm2]')


    dot.edge('interpolation_paste', 'param paste strength 28d')


    dot.edge('concrete homogenization', 'concrete strength 28d')
    dot.edge('aggregate content','volume computation')
    dot.node('geometry', 'cross section: width, height\n'
                         'beam span (lenght)', color=input)
    dot.edge('geometry', 'fem model')
    dot.edge('homogenized values', 'fem model')
    dot.node('homogenized values', 'concrete: density, nu, concrete thermal conductivity (kappa) \nconcrete heat capacity (C), heat release [J/m^3] \n concrete E 28 days [N/mm2]')
    dot.edge('concrete homogenization', 'homogenized values')
    dot.edge('concrete strength 28d', 'fem model')



    dot.node('calorimetry data','calorimetry data', color=input)
    dot.edge('calorimetry data','hydration identifcation')
    tex_node(dot,'hydration identifcation',r'parameter identification: $\phi_{\mathrm{hydr}}$', color=process)
    dot.edge('hydration identifcation', 'phi')


    tex_node(dot,'phi', r'$\phi_{\mathrm{hydr}}$')



    dot.node('strength data','concrete strength data', color=input)
    dot.edge('strength data','paste identifcation')
    #dot.node('strength identifcation','strength parameter\nidentification', color=process)


    tex_node(dot,'paste identifcation',r'parameter identification: $\phi_{\mathrm{paste}}$', color=process)



    dot.edge('paste identifcation', 'phi_paste')

    tex_node(dot,'phi_paste',r'$\phi_{\mathrm{paste}}$')

    dot.edge('phi_paste', 'interpolation_paste')




    dot.node('interpolation_paste','function to "interpolate" paste strength and E \nf(phi_paste,slag ratio)', color=process)




    dot.node('Youngs modulus function','Youngs modulus function', color=process)

    dot.edge('strength data','Youngs modulus function')
    dot.edge('Youngs modulus function', 'E')


    dot.node('E','approx concrete E data')
    dot.edge('E','paste identifcation')


    dot.edge('phi','interpolation')
    dot.edge('max doh','hydration identifcation')
    dot.node('interpolation','function to "interpolate" paramters \nf(phi,slag ratio)', color=process)
    dot.edge('interpolation','hydration parameters')
    dot.edge('interpolation','heat release binder')
    dot.node('hydration parameters', 'hydration parameters: \neta, B1, B2, T_ref, E_act')
    dot.node('heat release binder', 'heat release binder:\nQ_pot in J/g')
    dot.edge('hydration parameters','fem model')
    dot.edge('heat release binder','concrete homogenization')
    dot.edge('ratio_cemI_cemII','interpolation')


    dot.node('steel', 'steel: diameter [mm],\n n bars', color=input)
    dot.edge('steel', 'load bearing capacity')
    dot.edge('steel', 'co2 computation volume')

    dot.node('input beam', 'steel yield strenght [N/mm2]\n'
                           'concrete: cover [mm]', color=input)
    dot.edge('input beam', 'load bearing capacity')


    dot.node('fem input',' * initial and boundary temperature \n'
                         '* parameter ft(DoH): a_ft\n'
                         '* parameter fc(DoH): a_fc\n'
                         '* parameter E(DoH): alpha_t, alpha_0, a_E', color=input)
    dot.edge('fem input','fem model')

    dot.node('fem control','Control paramters:\n'
                           'dt, max_time', color=input)
    dot.edge('fem control','fem model')


    dot.node('legend_input_constant','Input (constant)', color=input)
    dot.node('legend_input_variable','Input (variable)', color=input, shape='rectangle')
    dot.node('legend_process','Process', color=process)
    dot.node('legend_parameter','Parameter')
    dot.node('legend_output', 'Output', color=kpi, shape='rectangle')


    dot.render(view=view, cleanup=True)


if __name__ == "__main__":
    paper_workflow_graph(file_name='test_output', view=True)
