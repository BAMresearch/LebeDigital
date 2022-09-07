import graphviz

process = 'blue'
input = 'green'
problem = 'red'
kpi = 'orange'


dot = graphviz.Digraph('LebeDigital_optimization_graph_ignored', comment='LebeDigital', format='pdf')


dot.node('cemI', 'cem I content \n[kg/m³]', color=input, shape='rectangle')
dot.edge('cemI','volume computation')
dot.edge('cemII','volume computation')
dot.edge('w','volume computation')
dot.edge('volume computation','ratio_cemI_cemII')

dot.node('cemII', 'cem II content \n[kg/m³]', color=input, shape='rectangle')


dot.node('plasticizer', 'plasticizer content [kg/m³]', color=input, shape='rectangle')


dot.node('aggregate content', 'aggregate content \n[kg/m³]', color=input, shape='rectangle')


dot.node('ratio_cemI_cemII', 'ratio of cem I/II [-]')



dot.node('w', 'water content\n w [kg/m³]', color=input, shape='rectangle')

dot.edge('volume computation','w/c')
dot.node('w/c', 'water cement ratio\n w/c [-]')





dot.edge('plasticizer', 'paste strength 28d')
dot.edge('plasticizer', 'E paste')



dot.node('volume computation', 'computation of volume contents\nand ratios', color=process)

dot.node('aggregate volume content', 'aggregate volume content []')
dot.edge('aggregate volume content', 'concrete homogenization')
dot.edge('volume computation','aggregate volume content')

dot.edge('ratio_cemI_cemII', 'E paste')

dot.edge('ratio_cemI_cemII', 'paste strength 28d')



dot.edge('w/c', 'E paste')
dot.edge('w/c', 'paste strength 28d')


dot.node('concrete E 28d', 'concrete E 28 days \n[N/mm2]')
dot.edge('concrete E 28d', 'load bearing capacity')

dot.node('concrete strength 28d', 'concrete strength \n28 days [N/mm2]')

dot.edge('ft(fc)', 'concrete tensile strength 28d')
dot.node('concrete tensile strength 28d', 'concrete tensile strength \n28 days [N/mm2]')
dot.node('ft(fc)', 'function ft(fc)', color=process)

dot.edge('concrete strength 28d', 'load bearing capacity')
dot.edge('concrete strength 28d', 'ft(fc)')

dot.node('cross section', 'cross section properties \n[mm]')
dot.edge('cross section', 'load bearing capacity')



dot.node('load bearing capacity', 'design code', color= process)
dot.edge('load bearing capacity', 'kpi load bearing capacity')
dot.node('kpi load bearing capacity', 'maximum traffic load \nin center 28d [kN]', color=kpi, shape='rectangle')

# new part from Erik


dot.node('demolding time', 'time of demolding [h]', color=kpi, shape='rectangle')
dot.node('fem model', 'fem model', color=process)
dot.edge('fem model','demolding time')


dot.node('concrete homogenization', 'concrete homogenization', color=process)
dot.edge('concrete tensile strength 28d','fem model')

dot.node('max temp', 'temperature constraint', color=input)
dot.edge('max temp','demolding time')


dot.node('ft paramter', 'parameter ft(DoH): a_ft', color=input)
dot.edge('ft paramter', 'fem model')

dot.node('fc paramter', 'parameter fc(DoH): a_fc', color=input)
dot.edge('fc paramter', 'fem model')


dot.node('E paramter', 'parameter E(DoH): \nalpha_t, alpha_0, a_E', color=input)
dot.edge('E paramter','fem model')


# just a test comment
dot.node('test', 'test', color=input)

dot.node('aggregate data', 'aggregate data: E, nu, rho, C, kappa', color=input)
dot.edge('aggregate data', 'concrete homogenization')


dot.node('paste data', 'paste data: nu, rho, C, kappa', color=input)
dot.edge('paste data', 'concrete homogenization')


dot.node('paste strength 28d', 'GP(paste strength after 28 days) \n[N/mm2]', color=process)
dot.node('param paste strength 28d', 'paste strength after 28 days \n[N/mm2]')

dot.edge('paste strength 28d', 'param paste strength 28d')
dot.edge('param paste strength 28d','concrete homogenization')
dot.edge('concrete homogenization', 'concrete strength 28d')
dot.edge('concrete homogenization', 'concrete E 28d')

dot.node('E paste', 'GP(paste Youngs modulus 28d) \n[N/mm2]', color=process)
dot.edge('E paste', 'param E paste')
dot.node('param E paste', 'paste Youngs modulus 28d \n[N/mm2]')
dot.edge('param E paste', 'concrete homogenization')

dot.edge('aggregate content','volume computation')

dot.node('geometry', 'part geometry', color=input, shape='rectangle')
dot.edge('geometry', 'fem model')
dot.edge('geometry', 'cross section')


dot.node('fem temperature','initial and boundary temperature', color=input, shape='rectangle')
dot.edge('fem temperature','fem model')
dot.edge('concrete density', 'fem model')
dot.node('concrete density', 'concrete density')
dot.edge('concrete density', 'fem model')
dot.node('concrete nu', 'concrete nu')
dot.node('concrete thermal', 'concrete thermal conductivity (kappa) \nconcrete heat capacity (C)')
dot.edge('concrete homogenization', 'concrete density')
dot.edge('concrete homogenization', 'concrete nu')
dot.edge('concrete homogenization', 'concrete thermal')

dot.edge('concrete thermal', 'fem model')

dot.edge('concrete nu', 'fem model')
dot.edge('concrete strength 28d', 'fem model')
dot.edge('concrete E 28d', 'fem model')

dot.node('phi', 'phi', color=input)
dot.edge('phi','interpolation')

dot.node('interpolation','function to "interpolate" paramters \nf(phi,cem ratio)', color=process)
dot.edge('interpolation','hydration parameters')

dot.node('hydration parameters', 'hydration parameters: \neta, B1, B2, Q_pot, T_ref, E_act',)
dot.edge('hydration parameters','fem model')
dot.node('time vector','time vector', color=input)
dot.edge('time vector','hydration model')

dot.node('reaction temperature','reaction temperature', color=input)
dot.edge('reaction temperature','hydration model')

dot.edge('hydration parameters','hydration model')
dot.node('hydration model','hydration model', color=process)
dot.edge('hydration model','hydration output')
dot.node('hydration output','data: heat over time', shape='rectangle')

dot.edge('ratio_cemI_cemII','interpolation')

dot.node('fkt_alpha_max','function: DoH_max(w/c)', color=process)
dot.edge('fkt_alpha_max','alpha_max')
dot.edge('w/c','fkt_alpha_max')
dot.node('alpha_max','max degree of hydration')
dot.edge('alpha_max', 'hydration model')
dot.edge('alpha_max', 'fem model')

dot.render(view=True)
