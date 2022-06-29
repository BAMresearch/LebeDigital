import graphviz
dot = graphviz.Digraph('LebeDigital', comment='LebeDigital ')


dot.node('cemII', 'cem II content c [kg/m³]')

dot.node('cemIII', 'cem III content c [kg/m³]')

dot.node('c', 'total cement content c [kg/m³]')

dot.edge('cemII', 'c')
dot.edge('cemIII', 'c')

dot.node('ratio_cem3_cem2', 'ratio of cem II/III\n r_II/III [-]')

dot.edge('cemII', 'ratio_cem3_cem2')
dot.edge('cemIII', 'ratio_cem3_cem2')

dot.node('w', 'water content\n w [kg/m³]')

dot.node('w/c', 'water cement ratio\n w/c [-]')

dot.edge('w', 'w/c')

dot.edge('c', 'w/c')

dot.render(directory='doctest-output', view=True)



#dot.edge('c', 'wc', constraint='false')
