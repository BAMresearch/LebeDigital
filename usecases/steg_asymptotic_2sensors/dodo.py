from pathlib import Path
from doit.task import clean_targets
import sys

DOIT_CONFIG = {"default_tasks": ["test"]}
PYTHON_EXE = sys.executable


def task_generate_virtual_metadata():
    """
    Create the metadata file for virtual experiments to be performed, store in yaml files.
    """
    script = Path(__file__).parents[0] / "generate_metadata_virtual_experiment.py"
    metadata_files = Path(__file__).parent.glob('*_meta.yaml')

    return {
        "actions": [f"{PYTHON_EXE} {script}"],
        "file_dep": [script],
        "verbosity": 2,  # show stdout
        "targets": [*metadata_files],
        'uptodate': [len(list(metadata_files)) > 0],
    }


def task_generate_virtual_samples():
    """
    Create virtual experiments to be stored in yaml files.
    """
    metadata_files = Path(__file__).parent.glob('*_meta.yaml')
    data_files = Path(__file__).parent.glob('*_data.yaml')

    script = Path(__file__).parents[0] / "virtual_experiment.py"

    return {
        "actions": [f"{PYTHON_EXE} {script}"],
        "file_dep": [script, *metadata_files],
        "verbosity": 2,  # show stdout
        "targets": [*data_files],
        "setup": ["generate_virtual_metadata"],
        "clean": [clean_targets]
    }


def task_optimize_asymptotic_model():
    """
    Optimize asymptotic model based on experimental data stored in yaml files.
    """
    script = Path(__file__).parents[0] / "test_optimize.py"
    script_dep_task = Path(__file__).parents[0] / "virtual_experiment.py"
    script_asymptotic_model = Path(__file__).parents[0] / "asymptotic_model.py"
    script_asymptotic_model_error = Path(__file__).parents[0] / "asymptotic_model_error.py"
    metadata_files = Path(__file__).parent.glob("*_meta.yaml")
    data_files = Path(__file__).parent.glob("*_data.yaml")

    return {
        "actions": [f"{PYTHON_EXE} {script}"],
        "file_dep": [script, script_dep_task, script_asymptotic_model, script_asymptotic_model_error,
                     *metadata_files, *data_files],
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
        "task_dep": ["optimize_asymptotic_model"]
    }
