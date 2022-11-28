import pathlib
from workflow_graph.paper_workflow_graph import paper_workflow_graph
from tex.macros.py_macros import py_macros
import yaml
from doit import get_var

DOIT_CONFIG = {
    "verbosity": 2,
}

# here is the possibility to implement different doit runs
# currently the default is set to 'mode=CI'
config = {"mode": get_var('mode', 'CI')}

# the mode or other variables can be used via a dictionary
# config['mode']

# general
ROOT = pathlib.Path(__file__).parent
figures_dir = 'figures'

# paper
paper_dir = 'tex'
paper_name = 'optimization_paper.tex'
paper_file = ROOT / paper_dir / paper_name


# macros
py_macros_file = ROOT / paper_dir / 'macros' / 'py_macros.tex'


# read macros yaml to define figure file names
with open(py_macros_file.with_suffix('.yaml')) as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)


# workflow graph
workflow_graph_dir = "workflow_graph"
workflow_graph_name = data['file_names']['workflowGraph']  # name of output pdf file as defined in macros yaml
workflow_graph_script = ROOT / workflow_graph_dir / 'paper_workflow_graph.py'
workflow_output_file = ROOT / figures_dir / workflow_graph_name


def task_build_graph():
    """build workflow graph"""
    return {
        "file_dep": [workflow_graph_script],
        "actions": [(paper_workflow_graph,[workflow_output_file.with_suffix(''),False])],
        "targets": [workflow_output_file],
        "clean": True,
    }


def task_build_tex_macros():
    """build tex macros"""
    return {
        "file_dep": [py_macros_file.with_suffix('.yaml')],
        "actions": [(py_macros, [py_macros_file.with_suffix('')])],
        "targets": [py_macros_file],
        "clean": True,
    }


def task_paper():
    """compile pdf from latex source"""
    return {
        "file_dep": [paper_file, workflow_output_file, py_macros_file],
        "actions": [f"tectonic {paper_file}"],
        "targets": [paper_file.with_suffix('.pdf')],
        "clean": True,
    }
