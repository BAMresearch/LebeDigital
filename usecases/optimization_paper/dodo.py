import pathlib
from workflow_graph.paper_workflow_graph import build_workflow_graph
from doit import create_after

# testing stuff
ROOT = pathlib.Path(__file__).parent

# directories
paper_dir = 'tex'
figures_dir = 'figures'

# paper
paper_name = 'optimization_paper'

# workflow graph
workflow_graph_script = "paper_workflow_graph.py"
workflow_graph_name = "optimization_paper_graph"
workflow_graph_dir = "workflow_graph"
workflow_output_name = figures_dir + '/' + workflow_graph_name
workflow_output_file = ROOT / figures_dir  / (workflow_graph_name + '.gz.pdf')



def task_build_graph():
    return {
        "file_dep": [workflow_graph_dir + '/' + workflow_graph_script],
        "actions": [(build_workflow_graph,[workflow_output_name,False])],
        "targets": [workflow_output_file],
        "clean": True,
    }


@create_after(executed='build_graph')
def task_paper():
    """compile pdf from latex source"""
    paper_tex = paper_dir + '/' + paper_name + '.tex'
    return {
        # TODO find out why this dependency throws an error
        "file_dep": [paper_tex, workflow_output_file],
        "actions": [f"tectonic {paper_tex}"],
        "targets": [paper_name + '.pdf'],
        "clean": True,
    }
