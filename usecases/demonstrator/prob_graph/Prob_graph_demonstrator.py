import graphviz

process = 'blue'
input = 'green'
problem = 'red'


dot = graphviz.Digraph('LebeDigital', comment='LebeDigital ')


dot.node('cemI', 'cem I content \n[kg/m³]', color=input)

dot.node('cemII', 'cem II content \n[kg/m³]', color=input)

dot.node('c', 'total cement content c \n[kg/m³]')

dot.node('plasticizer', 'plasticizer content [kg/m³]', color=input)


dot.node('aggregate content', 'aggregate content \n[kg/m³]', color=input)

dot.edge('cemI', 'c')
dot.edge('cemII', 'c')

dot.node('ratio_cemI_cemII', 'ratio of cem I/II [-]')
dot.edge('ratio_cemI_cemII','paste data')



dot.edge('cemI', 'ratio_cemI_cemII')
dot.edge('cemII', 'ratio_cemI_cemII')

dot.node('w', 'water content\n w [kg/m³]', color=input)

dot.node('w/c', 'water cement ratio\n w/c [-]', color=problem)

dot.edge('w', 'w/c')

dot.edge('c', 'w/c')

#dot.node('CO2_cemI', 'CO2 emission cemI [kg/m3]')
#dot.node('CO2_cemII', 'CO2 emission cemII [kg/m3]')
#dot.edge('CO2_cemI', 'CO2')
#dot.edge('CO2_cemII', 'CO2')

#dot.node('volume', 'volume of specimen [m3]')

#dot.edge('volume', 'CO2')

#dot.node('CO2', 'CO2 emission [kg]')

#dot.node('costs', 'cost [Euro]')


dot.edge('plasticizer', 'paste strength 28d')
dot.edge('plasticizer', 'E paste')


dot.node('paste/aggregate volume ratio', 'paste aggregate volume ratio []')

#dot.edge('aggregate content', '')
dot.edge('ratio_cemI_cemII', 'E paste')

dot.edge('ratio_cemI_cemII', 'paste strength 28d')



dot.edge('w/c', 'E paste')
dot.edge('w/c', 'paste strength 28d')


dot.node('concrete E 28d', 'concrete E 28 days \n[N/mm2]')
dot.edge('concrete E 28d', 'load bearing capacity')

dot.node('concrete strength 28d', 'concrete strength \n28 days [N/mm2]')
dot.node('concrete tensile strength 28d', 'concrete tensile strength \n28 days [N/mm2]')
dot.edge('concrete strength 28d', 'load bearing capacity')
dot.edge('concrete strength 28d', 'concrete tensile strength 28d')

dot.node('cross section', 'cross section properties \n[mm]')
dot.edge('cross section', 'load bearing capacity')



dot.node('load bearing capacity', 'maximum traffic load \nin center 28d [kN]')

# new part from Erik


dot.node('demolding time', 'time of demolding [h]')
dot.node('fem model', 'structural model', color=process)
dot.edge('fem model','demolding time')


dot.node('concrete homogenization', 'concrete homogenization', color=process)
dot.node('hydration model', 'hydration model', color=process)
dot.node('mechanics model', 'mechanics model', color=process)
dot.node('E(DoH)', 'function E(DoH)', color=process)
dot.node('fc(DoH)', 'function fc(DoH)', color=process)
dot.node('ft(DoH)', 'function ft(DoH)', color=process)
dot.edge('concrete tensile strength 28d','ft(DoH)')
dot.edge('ft(DoH)', 'mechanics model')

dot.node('max temp', 'temperature constraint', color=input)
dot.edge('max temp','demolding time')


dot.node('ft paramter', 'parameter ft(DoH): a_ft', color=problem)
dot.edge('ft paramter', 'ft(DoH)')

dot.node('cem I fc paramter', 'cem I parameter fc(DoH): a_fc', color=input)
dot.node('cem II fc paramter', 'cem II parameter fc(DoH): a_fc', color=input)
dot.edge('cem I fc paramter','mix fc fkt')
dot.edge('cem II fc paramter', 'mix fc fkt')
dot.node('mix fc fkt','interpolation or GP???',color=problem)
dot.edge('mix fc fkt', 'fc(DoH)')


dot.node('cem I E paramter', 'cem I parameter E(DoH): \nalpha_t, alpha_0, a_E', color=input)
dot.node('cem II E paramter', 'cem II parameter E(DoH): \nalpha_t, alpha_0, a_E', color=input)
dot.edge('cem I E paramter','mix E fkt')
dot.edge('cem II E paramter', 'mix E fkt')
dot.node('mix E fkt','interpolation or GP???',color=problem)
dot.edge('mix E fkt', 'E(DoH)')




dot.edge('hydration model', 'fem model')
dot.edge('mechanics model', 'fem model')
dot.edge('E(DoH)', 'mechanics model')
dot.edge('fc(DoH)', 'mechanics model')


dot.node('aggregate data', 'aggregate data: E, nu, rho, C, kappa', color=input)
dot.edge('aggregate data', 'concrete homogenization')

dot.node('cem I paste data', 'cem I paste data: nu, rho, C, kappa', color=input)
dot.edge('cem I paste data', 'paste data')
dot.node('cem II paste data', 'cem II paste data: nu, rho, C, kappa', color=input)
dot.edge('cem II paste data', 'paste data')

dot.node('paste data', 'paste data: nu, rho, C, kappa', color=process)
dot.edge('paste data', 'concrete homogenization')


dot.node('paste strength 28d', 'GP(paste strength after 28 days) \n[N/mm2]', color=process)
dot.edge('paste strength 28d', 'concrete homogenization')
dot.edge('concrete homogenization', 'concrete strength 28d')
dot.edge('concrete homogenization', 'concrete E 28d')

dot.node('E paste', 'GP(paste Youngs modulus 28d) \n[N/mm2]', color=process)
dot.edge('E paste', 'concrete homogenization')

dot.edge('aggregate content', 'paste/aggregate volume ratio')
dot.edge('paste/aggregate volume ratio', 'concrete homogenization')

dot.node('geometry', 'part geometry', color=input)
dot.edge('geometry', 'fem model')
dot.edge('geometry', 'cross section')


dot.node('concrete density', 'concrete density')
dot.node('concrete nu', 'concrete nu')
dot.node('concrete thermal', 'concrete thermal conductivity (kappa) \nconcrete heat capacity (C)')
dot.edge('concrete homogenization', 'concrete density')
dot.edge('concrete homogenization', 'concrete nu')
dot.edge('concrete homogenization', 'concrete thermal')


dot.edge('concrete density', 'hydration model')
dot.edge('concrete thermal', 'hydration model')


dot.edge('concrete nu', 'mechanics model')
dot.edge('concrete density', 'load')
dot.edge('concrete strength 28d', 'fc(DoH)')
dot.edge('concrete E 28d', 'E(DoH)')

dot.node('cem I hydration parameters', 'cem I hydration parameters: \neta, B1, B2, Q_pot, T_ref, E_act', color=input)
dot.edge('cem I hydration parameters','interpolation')
dot.node('cem II hydration parameters', 'cem II hydration parameters: \neta, B1, B2, Q_pot, T_ref, E_act', color=input)
dot.edge('cem II hydration parameters','interpolation')
dot.node('interpolation','interpolation or parallel models', color=problem)
dot.edge('interpolation','hydration model')

dot.edge('w/c', 'missing')
dot.node('missing', 'analytical function or GP???', color=problem)
dot.edge('missing','alpha max')
dot.node('alpha max', 'max degree of hydration')
dot.edge('alpha max','hydration model')

dot.node('load', 'load at demolding')

dot.edge('load','fem model')


dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
