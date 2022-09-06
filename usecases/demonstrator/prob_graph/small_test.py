import graphviz

process = 'blue'
input = 'green'
problem = 'red'
kpi = 'orange'


dot = graphviz.Digraph('LebeDigital', comment='LebeDigital',format='pdf')


dot.node('cemI', 'cem I content \n[kg/m³]', color=input, shape='rectangle')
dot.edge('cemI','volume computation')
dot.edge('cemII','volume computation')
dot.edge('w','volume computation')
dot.edge('volume computation','ratio_cemI_cemII')


dot.render(view=True)

#dot.render(view=True)
