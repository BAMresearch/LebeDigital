import graphviz

def paper_workflow_graph(file_name = 'test_output', view=False):
    print('test')

    file_name = str(file_name) # convert pathlib object to useful string

    process = 'blue'
    input = 'green'
    problem = 'red'
    kpi = 'orange'

    dot = graphviz.Digraph(file_name, comment='Optimization Paper', format='pdf')

    dot.node('epd products', 'environmental product decleration (epd) for cemI, cemII, aggregates in CO2/kg', color=input)
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
    dot.node('cemI', 'cem I, cem II, aggregate content \n[kg/m³]')
    dot.edge('volume computation', 'cemI')
    dot.edge('ratio_cemI_cemII','volume computation')
    dot.node('aggregate content', 'aggregate/(water+cement volume ratio) [-]', color=input, shape='rectangle')
    dot.node('ratio_cemI_cemII', 'volume ratio of cem I/II [-]', color=input, shape='rectangle')
    dot.edge('w/c','volume computation')
    dot.node('w/c', 'water cement mass ratio [-]\n plasticizer content [kg/m³]', color=input)
    dot.node('volume computation', 'computation of volume contents', color=process)
    dot.edge('aggregate content', 'concrete homogenization')
    dot.edge('ratio_cemI_cemII', 'E paste')
    dot.edge('ratio_cemI_cemII', 'paste strength 28d')
    dot.edge('w/c', 'E paste')
    dot.edge('w/c', 'paste strength 28d')
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
    dot.node('demolding time', 'time of demolding [h]', color=kpi, shape='rectangle')
    dot.node('fem model', 'fem model', color=process)
    dot.edge('fem model','demolding time')
    dot.node('concrete homogenization', 'concrete homogenization', color=process)
    dot.edge('concrete tensile strength 28d','fem model')
    dot.node('max temp', 'temperature constraint', color=input)
    dot.edge('max temp','demolding time')
    dot.node('aggregate data', 'aggregate data: E, nu, C, kappa', color=input)
    dot.node('aggregate rho', 'aggregate density', color=input)
    dot.node('cem rho', 'densities: cemI, cemII, water', color=input)
    dot.edge('cem rho','volume computation')
    dot.edge('volume computation','paste rho')
    dot.node('paste rho', 'paste density')
    dot.edge('paste rho',  'concrete homogenization')
    dot.edge('aggregate data', 'concrete homogenization')
    dot.edge('aggregate rho', 'volume computation')
    dot.edge('aggregate rho',  'concrete homogenization')
    dot.node('paste data', 'paste data: nu, C, kappa', color=input)
    dot.edge('paste data', 'concrete homogenization')
    dot.node('paste strength 28d', 'GP(paste strength after 28 days) \n[N/mm2]', color=process)
    dot.node('param paste strength 28d', 'paste strength after 28 days \n[N/mm2]')
    dot.edge('paste strength 28d', 'param paste strength 28d')
    dot.edge('param paste strength 28d','concrete homogenization')
    dot.edge('concrete homogenization', 'concrete strength 28d')
    dot.node('E paste', 'GP(paste Youngs modulus 28d) \n[N/mm2]', color=process)
    dot.edge('E paste', 'param E paste')
    dot.node('param E paste', 'paste Youngs modulus 28d \n[N/mm2]')
    dot.edge('param E paste', 'concrete homogenization')
    dot.edge('aggregate content','volume computation')
    dot.node('geometry', 'cross section: width, height\n'
                         'beam span (lenght)', color=input)
    dot.edge('geometry', 'fem model')
    dot.edge('homogenized values', 'fem model')
    dot.node('homogenized values', 'concrete: density, nu, concrete thermal conductivity (kappa) \nconcrete heat capacity (C), heat release [J/m^3] \n concrete E 28 days [N/mm2]')
    dot.edge('concrete homogenization', 'homogenized values')
    dot.edge('concrete strength 28d', 'fem model')
    dot.node('phi', 'phi', color=input)
    dot.edge('phi','interpolation')
    dot.node('interpolation','function to "interpolate" paramters \nf(phi,cem ratio)', color=process)
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

    dot.render(view=view)

if __name__ == "__main__":
    paper_workflow_graph(file_name='test_output', view=True)
