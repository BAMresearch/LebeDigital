import pathlib
from workflow_graph.paper_workflow_graph import paper_workflow_graph
from tex.macros.py_macros import py_macros
import yaml
from doit import get_var
from doit.action import CmdAction

DOIT_CONFIG = {
    "verbosity": 2,
}

# here is the possibility to implement different doit runs
# currently the default is set to 'mode=default'
config = {"mode": get_var('mode', 'default')}

# the mode or other variables can be used via a dictionary
# config['mode']

# general
ROOT = pathlib.Path(__file__).parent
figures_dir = 'figures'

# initialize a list of file dependencies for the paper
PAPER_PLOTS = []

# function to return target path and append paper dependency list
def paper_plot_target(name):
    target = ROOT / figures_dir / name
    PAPER_PLOTS.append(target)
    return target


# paper
paper_dir = 'tex'
paper_name = 'optimization_paper.tex'
paper_file = ROOT / paper_dir / paper_name


# macros
py_macros_file = ROOT / paper_dir / 'macros' / 'py_macros.tex'
tex_macros_file = ROOT / paper_dir / 'macros' / 'tex_macros.tex'
TEX_MACROS = [py_macros_file, tex_macros_file]


# read macros yaml to define figure file names
with open(py_macros_file.with_suffix('.yaml')) as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)


# workflow graph
workflow_graph_dir = "workflow_graph"
workflow_graph_name = data['file_names']['workflowGraph']  # name of output pdf file as defined in macros yaml
workflow_graph_script = ROOT / workflow_graph_dir / 'paper_workflow_graph.py'
workflow_output_file = ROOT / figures_dir / workflow_graph_name


def task_build_snakemake_dag():
    """build snakemake optimization workflow graph"""
    output_file_name = data['file_names']['snakemakeGraph']  # name of output pdf file as defined in macros yaml
    snakemake_dir = 'optimization_workflow'
    snakefile = ROOT / snakemake_dir / 'Snakefile'

    target = paper_plot_target(output_file_name)

    return {
        "file_dep": [snakefile],
        "actions": [CmdAction(f'cd {snakemake_dir} && snakemake --forceall --dag | dot -Tpdf > {target}')],
        "targets": [target],
        "clean": True,
    }


def task_build_graph():
    """build workflow graph"""
    target = paper_plot_target(workflow_graph_name)

    return {
        "file_dep": [workflow_graph_script],
        "actions": [(paper_workflow_graph,[target.with_suffix('').with_suffix(''),False])],
        "targets": [target],
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
        "file_dep": [paper_file] + TEX_MACROS + PAPER_PLOTS,
        "actions": [f"tectonic {paper_file}"],
        "targets": [paper_file.with_suffix('.pdf')],
        "clean": True,
    }
