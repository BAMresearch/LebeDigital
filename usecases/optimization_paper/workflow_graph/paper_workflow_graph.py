import graphviz

def paper_workflow_graph(file_name = 'test_output', view=False):

    file_name = str(file_name) # convert pathlib object to useful string

    process = 'blue'
    input = 'green'
    problem = 'red'
    kpi = 'orange'

    dot = graphviz.Digraph(file_name, comment='Optimization Paper', format='pdf')

    dot.node('epd products', 'environmental product decleration (epd) for cement, slag, aggregates in CO2/kg', color=input)
    dot.edge('epd products','co2 computation')
    dot.node('co2 computation', 'computation of global warming potential (GWP) per m³', color=process)
    dot.edge('co2 computation','gwp')
    dot.node('gwp', 'GWP per m³ [kg CO2,eq/m³]')
    dot.edge('gwp','co2 computation volume')
    dot.node('co2 computation volume', 'GWP per part', color=process)
    dot.edge('geometry','CO2 computation volume')
    dot.edge('co2 computation volume','gwp total')
    dot.edge('cemI','co2 computation')
    dot.node('gwp total', 'GWP per part [kg CO2,eq]', color=kpi, shape='rectangle')
    dot.node('cemI', 'cement, slag, aggregate content \n[kg/m³]')
    dot.edge('volume computation', 'cemI')
    dot.edge('ratio_cemI_cemII','volume computation')
    dot.node('aggregate content', 'aggregate/(water+cement volume ratio) [-]', color=input, shape='rectangle')
    dot.node('ratio_cemI_cemII', 'slag to concrete ratio', color=input, shape='rectangle')
    dot.edge('w/c','volume computation')
    dot.node('w/c', 'water cement mass ratio [-]\n plasticizer content [kg/m³]', color=input)
    dot.node('volume computation', 'computation of volume contents', color=process)
    dot.edge('aggregate content', 'concrete homogenization')
    dot.edge('ratio_cemI_cemII',  'interpolation_E')


    dot.edge( 'phi_E_paste','interpolation_E')
    dot.edge( 'interpolation_E','param paste E 28d')

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
    dot.node('param paste strength 28d', 'paste strength after 28 days [N/mm2]')
    dot.node('param paste E 28d', 'paste Youngs modulus after 28 days [N/mm2]')
    dot.edge('param paste E 28d', 'concrete homogenization')


    dot.edge('interpolation_strength', 'param paste strength 28d')


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
    dot.node('hydration identifcation','hydration parameter\nidentification', color=process)
    dot.edge('hydration identifcation', 'phi')
    dot.node('phi', 'phi_hydration')



    dot.node('strength data','concrete strength data', color=input)
    dot.edge('strength data','strength identifcation')
    dot.node('strength identifcation','strength parameter\nidentification', color=process)
    dot.edge('strength identifcation', 'phi_strength')
    dot.node('phi_strength', 'phi_strength')



    dot.edge('ratio_cemI_cemII','interpolation_strength')
    dot.edge('phi_strength','interpolation_strength')
    dot.node('interpolation_strength','function to "interpolate" paste strength \nf(phi_strength,slag ratio)', color=process)

    dot.node('interpolation_E','function to "interpolate" paste E \nf(phi_strength,slag ratio)', color=process)






    dot.node('Youngs modulus function','Youngs modulus function', color=process)



    dot.node('E identifcation','E_paste parameter\nidentification', color=process)
    dot.edge('E identifcation', 'phi_E_paste')
    dot.node('phi_E_paste', 'phi_E_paste')


    dot.edge('strength data','Youngs modulus function')
    dot.edge('Youngs modulus function', 'E')


    dot.node('E','approx concrete E data')
    dot.edge('E','E identifcation')


    dot.edge('phi','interpolation')
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
                         '* max degree of hydration \n'
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


    dot.render(view=view)

if __name__ == "__main__":
    paper_workflow_graph(file_name='test_output', view=True)
