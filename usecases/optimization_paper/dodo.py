import collections.abc
import pathlib

import yaml
from bam_figures.beam_design_plot import beam_design_plot
from bam_figures.create_heat_release_figure import create_heat_release_figure
from bam_figures.create_homogenization_figure import create_homogenization_figure
from bam_figures.create_mechanics_evolution_figure import create_mechanics_evolution_figure
from bam_figures.design_approach import design_approach_graph
from bam_figures.paper_workflow_graph import paper_workflow_graph
from doit import get_var
from doit.action import CmdAction
from tex.macros.py_macros import input_optimization_macros, py_macros


def update(d, u):
    """Merging nested dictionaries without deleting data
    from https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    Required to merge the two py marco yaml files
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


# function to return target path and append paper dependency list
def paper_plot_target(name):
    target = ROOT / FIGURES_DIR / name
    PAPER_PLOTS.append(target)
    return target


# config doit
DOIT_CONFIG = {"verbosity": 2}

# here is the possibility to implement different doit runs
# currently the default is set to 'mode=default'
# use mode=final to create the final paper with expensive figures etc.
# the mode or other variables can be used via a dictionary
# config['mode']
config = {"mode": get_var("mode", "default")}

# defining relevant paths and directories
ROOT = pathlib.Path(__file__).parent
FIGURES_DIR = "figures"  # target folder for all figures
BAM_PLOT_DIR = "bam_figures"  # origin folder for all BAM figure scripts
PAPER_DIR = pathlib.Path("tex")  # folder for the paper code
OPTIMIZATION_DIR = "optimization_workflow"  # folder for the paper code
MACROS_DIR = PAPER_DIR / "macros"  # folder for the paper code
path_to_optimization_workflows = ROOT / OPTIMIZATION_DIR

# initialize a list of file dependencies for the paper
PAPER_PLOTS = []

# macros
py_macros_file_BAM = ROOT / MACROS_DIR / "py_macros_BAM.tex"
py_macros_file_TUM = ROOT / MACROS_DIR / "py_macros_TUM.tex"
tex_macros_file_BAM = ROOT / MACROS_DIR / "tex_macros_BAM.tex"
tex_macros_file_TUM = ROOT / MACROS_DIR / "tex_macros_TUM.tex"
py_macros_optimization_file = ROOT / MACROS_DIR / "py_macros_optimization_input.tex"

TEX_MACROS = [
    py_macros_file_BAM,
    py_macros_file_TUM,
    tex_macros_file_BAM,
    tex_macros_file_TUM,
    py_macros_optimization_file,
]

# read macros yaml to define figure file names
with open(py_macros_file_TUM.with_suffix(".yaml")) as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)
with open(py_macros_file_BAM.with_suffix(".yaml")) as f:
    data = update(data, yaml.load(f, Loader=yaml.SafeLoader))


def task_build_graph():
    """build workflow graph"""
    # workflow graph
    workflow_graph_name = data["file_names"]["workflowGraph"]  # name of output pdf file as defined in macros yaml
    workflow_graph_script = ROOT / BAM_PLOT_DIR / "paper_workflow_graph.py"

    target = paper_plot_target(workflow_graph_name)

    return {
        "file_dep": [workflow_graph_script],
        "actions": [(paper_workflow_graph, [target.with_suffix("").with_suffix(""), False])],
        "targets": [target],
        "clean": True,
    }


def task_build_snakemake_dag():
    """build snakemake optimization workflow graph"""
    output_file_name = data["file_names"]["snakemakeGraph"]  # name of output pdf file as defined in macros yaml
    snakefile = ROOT / OPTIMIZATION_DIR / "Snakefile"

    target = paper_plot_target(output_file_name)

    return {
        "file_dep": [snakefile],
        "actions": [
            CmdAction(f"cd {OPTIMIZATION_DIR} && snakemake --forceall --dag 2> log_dag.txt | dot -Tpdf > {target}")
        ],
        "targets": [target],
        "clean": True,
    }


def task_build_design_figures():
    """build design approach graph"""

    design_apporach_plot_script = ROOT / BAM_PLOT_DIR / "design_approach.py"

    output_file_name_1 = data["file_names"]["designStandard"]  # name of output pdf file as defined in macros yaml
    output_file_name_2 = data["file_names"]["designProposed"]  # name of output pdf file as defined in macros yaml

    target1 = paper_plot_target(output_file_name_1)
    target2 = paper_plot_target(output_file_name_2)

    return {
        "file_dep": [design_apporach_plot_script],
        "actions": [
            (
                design_approach_graph,
                [target1.with_suffix("").with_suffix(""), target2.with_suffix("").with_suffix("")],
            )
        ],
        "targets": [target1, target2],
        "clean": True,
    }


def task_build_homogenization_figure():
    """build homogenization figure"""
    # homogenization figure
    homogenization_plot_name = data["file_names"]["homogenizationPlot"]  # name of pdf file defined in macros yaml
    homogenization_plot_script = ROOT / BAM_PLOT_DIR / "create_homogenization_figure.py"
    homogenization_plot_output_file = ROOT / FIGURES_DIR / homogenization_plot_name

    target = paper_plot_target(homogenization_plot_name)

    return {
        "file_dep": [homogenization_plot_script, py_macros_file_BAM.with_suffix(".yaml")],
        "actions": [
            (
                create_homogenization_figure,
                [data["homogenization_example_parameters"], homogenization_plot_output_file],
            )
        ],
        "targets": [target],
        "clean": True,
    }


def task_build_heat_release_figure():
    """build homogenization figure"""
    # homogenization figure
    heat_release_plot_name = data["file_names"]["heatReleasePlot"]  # name of pdf file defined in macros yaml
    heat_release_plot_script = ROOT / BAM_PLOT_DIR / "create_heat_release_figure.py"
    heat_release_plot_output_file = ROOT / FIGURES_DIR / heat_release_plot_name

    target = paper_plot_target(heat_release_plot_name)

    return {
        "file_dep": [heat_release_plot_script, py_macros_file_BAM.with_suffix(".yaml")],
        "actions": [
            (
                create_heat_release_figure,
                [data["heat_release_example_parameters"], heat_release_plot_output_file],
            )
        ],
        "targets": [target],
        "clean": True,
    }


def task_build_evolution_figure():
    """build evolution figure"""
    # evolution figure
    evolution_plot_name = data["file_names"]["evolutionPlot"]  # name of pdf file defined in macros yaml
    evolution_plot_script = ROOT / BAM_PLOT_DIR / "create_mechanics_evolution_figure.py"
    evolution_plot_output_file = ROOT / FIGURES_DIR / evolution_plot_name

    target = paper_plot_target(evolution_plot_name)

    return {
        "file_dep": [evolution_plot_script, py_macros_file_BAM.with_suffix(".yaml")],
        "actions": [
            (
                create_mechanics_evolution_figure,
                [data["evolution_example_parameters"], evolution_plot_output_file],
            )
        ],
        "targets": [target],
        "clean": True,
    }


def task_build_beam_design_figure():
    """build beam design figure"""
    beam_design_plot_name = data["file_names"]["beamDesignPlot"]  # name of output pdf file as defined in macros yaml
    beam_design_plot_script = ROOT / BAM_PLOT_DIR / "beam_design_plot.py"
    beam_design_plot_output_file = ROOT / FIGURES_DIR / beam_design_plot_name

    target = paper_plot_target(beam_design_plot_name)

    if config["mode"] == "default":
        n = 10  # number of points for figure, 150 should be good for final paper, but takes too long for development
    elif config["mode"] == "final":
        n = 150
    else:
        raise ValueError("Unknown mode")

    return {
        "file_dep": [beam_design_plot_script, py_macros_file_BAM.with_suffix(".yaml")],
        "actions": [(beam_design_plot, [data["beam_design_example_parameters"], n, beam_design_plot_output_file])],
        "targets": [target],
        "clean": True,
    }


def task_build_tex_macros():
    """build tex macros"""
    return {
        "file_dep": [py_macros_file_BAM.with_suffix(".yaml"), py_macros_file_TUM.with_suffix(".yaml")],
        "actions": [(py_macros, [[py_macros_file_BAM.with_suffix(""), py_macros_file_TUM.with_suffix("")]])],
        "targets": [py_macros_file_BAM, py_macros_file_TUM],
        "clean": True,
    }


def task_build_optimization_tex_macros():
    """build optimization input and results tex macros"""
    script = ROOT / MACROS_DIR / "py_macros.py"

    return {
        "file_dep": [script],
        "actions": [(input_optimization_macros, [path_to_optimization_workflows, py_macros_optimization_file])],
        "targets": [py_macros_optimization_file],
        "clean": True,
    }


def task_paper():
    """compile pdf from latex source"""
    paper_name = "optimization_paper.tex"
    paper_file = ROOT / PAPER_DIR / paper_name

    return {
        "file_dep": [paper_file] + TEX_MACROS + PAPER_PLOTS,
        "actions": [f"tectonic {paper_file} -c minimal 2> log_tectonic.txt"],
        "targets": [paper_file.with_suffix(".pdf")],
        "clean": True,
    }
