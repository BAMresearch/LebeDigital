import graphviz
dot = graphviz.Digraph('LebeDigital', comment='LebeDigital ')


dot.node('cemI', 'cem I content c [kg/m³]')

dot.node('cemII', 'cem II content c [kg/m³]')

dot.node('c', 'total cement content c [kg/m³]')

dot.node('plasticizer', 'plasticizer content [kg/m³]')


dot.node('aggregate content', 'aggregate content [kg/m³]')

dot.edge('cemI', 'c')
dot.edge('cemII', 'c')

dot.node('ratio_cemI_cemII', 'ratio of cem I/II\n r_I/II [-]')

dot.edge('cemI', 'ratio_cemI_cemII')
dot.edge('cemII', 'ratio_cemI_cemII')

dot.node('w', 'water content\n w [kg/m³]')

dot.node('w/c', 'water cement ratio\n w/c [-]')

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


dot.node('concrete E 28d', 'concrete E 28 days [N/mm2]')
dot.edge('concrete E 28d', 'load bearing capacity')

dot.node('concrete strength 28d', 'concrete strength 28 days [N/mm2]')
dot.node('concrete tensile strength 28d', 'concrete tensile strength 28 days [N/mm2]')
dot.edge('concrete strength 28d', 'load bearing capacity')
dot.edge('concrete strength 28d', 'concrete tensile strength 28d')

dot.node('cross section', 'cross section properties [mm]')
dot.edge('cross section', 'load bearing capacity')



dot.node('load bearing capacity', 'maximum traffic load in center 28d [kN]')

# new part from Erik


dot.node('demolding time', 'time of demolding [h]')
dot.node('fem model', 'structural model')
dot.edge('fem model','demolding time')


dot.node('concrete homogenization', 'concrete homogenization')
dot.node('hydration model', 'hydration model')
dot.node('mechanics model', 'mechanics model')
dot.node('E(DoH)', 'function E(DoH)')
dot.node('fc(DoH)', 'function fc(DoH)')
dot.node('ft(DoH)', 'function ft(DoH)')
dot.edge('concrete tensile strength 28d','ft(DoH)')
dot.edge('ft(DoH)', 'mechanics model')

dot.node('max temp', 'temperature constraint')
dot.edge('max temp','demolding time')


dot.node('ft paramter', 'parameter ft(DoH): a_ft')
dot.edge('ft paramter', 'ft(DoH)')
dot.node('fc paramter', 'parameter fc(DoH): a_fc')
dot.edge('fc paramter', 'fc(DoH)')
dot.node('E paramter', 'parameter E(DoH): alpha_t, alpha_0, a_E')
dot.edge('E paramter', 'E(DoH)')

dot.edge('hydration model', 'fem model')
dot.edge('mechanics model', 'fem model')
dot.edge('E(DoH)', 'mechanics model')
dot.edge('fc(DoH)', 'mechanics model')


dot.node('aggregate data', 'aggregate data: E, nu, rho, C, kappa')
dot.edge('aggregate data', 'concrete homogenization')
dot.node('paste data', 'paste data: nu, rho, C, kappa')
dot.edge('paste data', 'concrete homogenization')


dot.node('Phi GP strength paste', 'Phi GP strength paste')
dot.edge('Phi GP strength paste','paste strength 28d')


dot.node('paste strength 28d', 'GP(paste strength after 28 days) [N/mm2]')
dot.edge('paste strength 28d', 'concrete homogenization')
dot.edge('concrete homogenization', 'concrete strength 28d')
dot.edge('concrete homogenization', 'concrete E 28d')

dot.node('E paste', 'GP(paste Youngs modulus 28d) [N/mm2]')
dot.edge('E paste', 'concrete homogenization')

dot.edge('aggregate content', 'paste/aggregate volume ratio')
dot.edge('paste/aggregate volume ratio', 'concrete homogenization')

dot.node('geometry', 'part geometry')
dot.edge('geometry', 'fem model')
dot.edge('geometry', 'cross section')


dot.node('concrete density', 'concrete density')
dot.node('concrete nu', 'concrete nu')
dot.node('concrete thermal', 'concrete thermal conductivity (kappa) / heat capacity (C)')
dot.edge('concrete homogenization', 'concrete density')
dot.edge('concrete homogenization', 'concrete nu')
dot.edge('concrete homogenization', 'concrete thermal')


dot.edge('concrete density', 'hydration model')
dot.edge('concrete thermal', 'hydration model')


dot.edge('concrete nu', 'mechanics model')
dot.edge('concrete density', 'load')
dot.edge('concrete strength 28d', 'fc(DoH)')
dot.edge('concrete E 28d', 'E(DoH)')

dot.node('hydration parameters', 'hydration parameters: eta, B1, B2, Q_pot, T_ref, E_act')
dot.edge('hydration parameters','hydration model')

dot.edge('w/c', 'missing')
dot.node('missing', '???????')
dot.edge('missing','alpha max')
dot.node('alpha max', 'max degree of hydration')
dot.edge('alpha max','hydration model')

dot.node('load', 'load at demolding')

dot.edge('load','fem model')


dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
