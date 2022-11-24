import pathlib
from workflow_graph.paper_workflow_graph import paper_workflow_graph
from doit import create_after



# general
ROOT = pathlib.Path(__file__).parent
figures_dir = 'figures'

# paper
paper_dir = 'tex'
paper_name = 'optimization_paper.tex'
paper_file = ROOT / paper_dir / paper_name

# workflow graph
workflow_graph_dir = "workflow_graph"
workflow_graph_name = "paper_workflow_graph"
workflow_graph_script = ROOT / workflow_graph_dir / 'paper_workflow_graph.py'
workflow_output_file = ROOT / figures_dir / workflow_graph_name

def task_build_graph():
    """build workflow graph"""
    return {
        "file_dep": [workflow_graph_script],
        "actions": [(paper_workflow_graph,[workflow_output_file,False])],
        "targets": [workflow_output_file.with_suffix('.gv.pdf')],
        "clean": True,
    }

def task_paper():
    """compile pdf from latex source"""
    return {
        # TODO find out why this dependency throws an error
        "file_dep": [paper_file, workflow_output_file.with_suffix('.gv.pdf')],
        "actions": [f"tectonic {paper_file}"],
        "targets": [paper_file.with_suffix('.pdf')],
        "clean": True,
    }
