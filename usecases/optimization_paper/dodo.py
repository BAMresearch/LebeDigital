import pathlib
from workflow_graph.paper_workflow_graph import paper_workflow_graph
from bam_figures.create_homogenization_figure import create_homogenization_figure
from bam_figures.beam_design_plot import beam_design_plot
from tex.macros.py_macros import py_macros
import yaml
from doit import get_var
from doit.action import CmdAction
import collections.abc


def update(d, u):
    """ Merging nested dictionaries without deleting data
        from https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
        Required to merge the two py marco yaml files
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


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
py_macros_file_BAM = ROOT / paper_dir / 'macros' / 'py_macros_BAM.tex'
py_macros_file_TUM = ROOT / paper_dir / 'macros' / 'py_macros_TUM.tex'
tex_macros_file_BAM = ROOT / paper_dir / 'macros' / 'tex_macros_BAM.tex'
tex_macros_file_TUM = ROOT / paper_dir / 'macros' / 'tex_macros_TUM.tex'
TEX_MACROS = [py_macros_file_BAM, py_macros_file_TUM, tex_macros_file_BAM, tex_macros_file_TUM]






# read macros yaml to define figure file names
with open(py_macros_file_TUM.with_suffix('.yaml')) as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)
with open(py_macros_file_BAM.with_suffix('.yaml')) as f:
    data = update(data, yaml.load(f, Loader=yaml.SafeLoader))


# workflow graph
workflow_graph_dir = "workflow_graph"
workflow_graph_name = data['file_names']['workflowGraph']  # name of output pdf file as defined in macros yaml
workflow_graph_script = ROOT / workflow_graph_dir / 'paper_workflow_graph.py'
workflow_output_file = ROOT / figures_dir / workflow_graph_name

def task_build_graph():
    """build workflow graph"""
    target = paper_plot_target(workflow_graph_name)

    return {
        "file_dep": [workflow_graph_script],
        "actions": [(paper_workflow_graph,[target.with_suffix('').with_suffix(''),False])],
        "targets": [target],
        "clean": True,
    }


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


# homogenization figure
bam_plot_dir = "bam_figures"
homogenization_plot_name = data['file_names']['homogenizationPlot']  # name of output pdf file as defined in macros yaml
homogenization_plot_script = ROOT / bam_plot_dir / 'create_homogenization_figure.py'
homogenization_plot_output_file = ROOT / figures_dir / homogenization_plot_name

def task_build_homogenization_figure():
    """build homogenization figure"""
    target = paper_plot_target(homogenization_plot_name)

    return {
        "file_dep": [homogenization_plot_script],
        "actions": [(create_homogenization_figure,[data['homogenization_example_parameters'], homogenization_plot_output_file])],
        "targets": [target],
        "clean": True,
    }

# homogenization figure
beam_design_plot_name = data['file_names']['beamDesignPlot']  # name of output pdf file as defined in macros yaml
beam_design_plot_script = ROOT / bam_plot_dir / 'beam_design_plot.py'
beam_design_plot_output_file = ROOT / figures_dir / beam_design_plot_name

def task_build_beam_design_figure():
    """build beam design figure"""
    target = paper_plot_target(beam_design_plot_name)
    n = 10  # number of points for figure, 150 should be good for final paper, but takes too long for development

    return {
        "file_dep": [beam_design_plot_script],
        "actions": [(beam_design_plot,[data['beam_design_example_parameters'], n, beam_design_plot_output_file])],
        "targets": [target],
        "clean": True,
    }

def task_build_tex_macros():
    """build tex macros"""
    return {
        "file_dep": [py_macros_file_BAM.with_suffix('.yaml'),py_macros_file_TUM.with_suffix('.yaml')],
        "actions": [(py_macros, [[py_macros_file_BAM.with_suffix(''),py_macros_file_TUM.with_suffix('')]])],
        "targets": [py_macros_file_BAM,py_macros_file_TUM],
        "clean": True,
    }


def task_paper():
    """compile pdf from latex source"""
    return {
        "file_dep": [paper_file] + TEX_MACROS + PAPER_PLOTS,
        "actions": [f"tectonic {paper_file}"], # AA: I need to add./tectonic
        "targets": [paper_file.with_suffix('.pdf')],
        "clean": True,
    }

