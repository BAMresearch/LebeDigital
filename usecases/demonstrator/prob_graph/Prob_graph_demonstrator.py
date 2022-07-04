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

dot.node('Phi GP strength paste', 'Phi GP strength paste')
dot.edge('Phi GP strength paste','paste strength 28d')


dot.node('paste strength 28d', 'GP(paste strength after 28 days) [N/mm2]')
dot.edge('paste strength 28d', 'concrete strength, E 28d')

dot.edge('aggregate content', 'paste/aggregate volume ratio')

dot.edge('plasticizer', 'paste strength 28d')
dot.edge('plasticizer', 'E paste')


dot.node('paste/aggregate volume ratio', 'paste aggregate volume ratio []')
dot.edge('paste/aggregate volume ratio', 'concrete strength, E 28d')

#dot.edge('aggregate content', '')
dot.edge('ratio_cemI_cemII', 'E paste')

dot.edge('ratio_cemI_cemII', 'paste strength 28d')


dot.node('E paste', 'GP(paste Youngs modulus 28d) [N/mm2]')
dot.edge('E paste', 'concrete strength, E 28d')

dot.edge('w/c', 'E paste')
dot.edge('w/c', 'paste strength 28d')


dot.node('nu paste, E,nu aggregates', 'paste nu and aggregate '
                                       'E/nu[-]')
dot.edge('nu paste, E,nu aggregates', 'concrete strength, E 28d')


dot.node('concrete strength, E 28d', 'concrete strength, E after 28 days ['
                                     'N/mm2]')
dot.edge('concrete strength, E 28d', 'load bearing capacity')


dot.node('cross section', 'cross section properties [mm]')
dot.edge('cross section', 'load bearing capacity')



dot.node('load bearing capacity', 'maximum traffic load in center 28d [kN]')


dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
