from pathlib import Path
from doit.task import clean_targets
import sys

DOIT_CONFIG = {"default_tasks": ["test"]}
PYTHON_EXE = sys.executable


def task_generate_virtual_samples():
    """
    Create virtual experiments to be stored in yaml files.
    """
    script = Path(__file__).parents[0] / "virtual_experiment.py"
    result_file_linear = Path(__file__).parents[0] / "virtual_linear_experiment_model.yaml"
    result_file_quadratic = Path(__file__).parents[0] / "virtual_quadratic_experiment_model.yaml"
    p = (Path(__file__).parents[0]).glob("virtual_*_experiment_data_*.yaml")
    result_files_data = [x for x in p if x.is_file()]

    return {
        "actions": [f"{PYTHON_EXE} {script}"],
        "file_dep": [script],
        "verbosity": 2,  # show stdout
        "targets": [result_file_linear, result_file_quadratic] + result_files_data,
        "clean": [clean_targets]
    }


def task_optimize_linear_model():
    """
    Optimize linear model based on experimental data stored in yaml files.
    """
    script = Path(__file__).parents[0] / "test_optimize.py"
    script_dep_task = Path(__file__).parents[0] / "virtual_experiment.py"
    script_linear_model = Path(__file__).parents[0] / "linear_model.py"
    script_linear_model_error = Path(__file__).parents[0] / "linear_model_error.py"

    return {
        "actions": [f"{PYTHON_EXE} {script}"],
        "file_dep": [script, script_dep_task, script_linear_model, script_linear_model_error],
        "setup": ["generate_virtual_samples"],
        "verbosity": 2,  # show stdout
    }


def task_flake():
    return {
        "actions": ["flake8 *.py --count --exit-zero --select=E9,F63,F7,F82 --show-source --statistics"],
        "file_dep": ["dodo.py"],
        "verbosity": 2,  # show stdout
    }


def task_test():
    """ group all tests """
    # specifying each test as a separate tasks allows to include all the dependencies separately
    return {
        "actions": None,
        "task_dep": ["optimize_linear_model"]
    }
