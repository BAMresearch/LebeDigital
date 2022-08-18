import graphviz

process = 'blue'
input = 'green'
problem = 'red'
kpi = 'orange'


dot = graphviz.Digraph('LebeDigital', comment='LebeDigital ')


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

dot.node('geometry', 'part geometry', color=input)
dot.edge('geometry', 'fem model')
dot.edge('geometry', 'cross section')


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

dot.node('cem I hydration parameters', 'cem I hydration parameters: \neta, B1, B2, Q_pot, T_ref, E_act', color=input)
dot.edge('cem I hydration parameters','interpolation')
dot.node('cem II hydration parameters', 'cem II hydration parameters: \neta, B1, B2, Q_pot, T_ref, E_act', color=input)
dot.edge('cem II hydration parameters','interpolation')
dot.node('interpolation','TODO:\ninterpolation or parallel models', color=problem)
dot.edge('interpolation','fem model')
dot.edge('ratio_cemI_cemII','interpolation')

dot.edge('w/c', 'fem model')



dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
