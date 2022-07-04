import graphviz
dot = graphviz.Digraph('LebeDigital', comment='LebeDigital '2
)


dot.node('cemI', 'cem I content c [kg/m³]')

dot.node('cemII', 'cem II content c [kg/m³]')

dot.node('c', 'total cement content c [kg/m³]')

dot.edge('cemI', 'c')
dot.edge('cemII', 'c')

dot.node('ratio_cemI_cemII', 'ratio of cem I/II\n r_I/II [-]')

dot.edge('cemI', '@')
dot.edge('cemII', 'ratio_cemI_cemII')

dot.node('w', 'water content\n w [kg/m³]')

dot.node('w/c', 'water cement ratio\n w/c [-]')

dot.edge('w', 'w/c')

dot.edge('c', 'w/c')

dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
