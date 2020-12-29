from pathlib import Path
from doit.task import clean_targets
import yaml

def task_generate_virtual_samples():
    """
    Create virtual experiments to be stored in yaml files.
    """
    script = Path(__file__).parents[0] / "virtual_experiment.py"
    result_file_linear = Path(__file__).parents[0] / "virtual_linear_experiment_model.yaml"
    result_file_quadratic = Path(__file__).parents[0] / "virtual_quadratic_experiment_model.yaml"
    p = (Path(__file__).parents[0]).glob("virtual_*_experiment_data_*.yaml")
    result_files_data = [x for x in p if x.is_file()]

    yield {
            "basename" : "TASK: generate virtual samples",
            "actions": [f"python {script}"],
            "file_dep" : [script],
            "verbosity": 2, # show stdout
            "targets": [result_file_linear, result_file_quadratic] + result_files_data,
            "clean": [clean_targets]
            }

def task_optimize_linear_model():
    """
    Optimize linear model based on experimental data stored in yaml files.
    """
    script = Path(__file__).parents[0] / "optimize.py"
    script_linear_model = Path(__file__).parents[0] / "linear_model.py"
    script_linear_model_error = Path(__file__).parents[0] / "linear_model_error.py"
    exp_linear = Path(__file__).parents[0] / "virtual_linear_experiment_model.yaml"
    exp_quadratic = Path(__file__).parents[0] / "virtual_quadratic_experiment_model.yaml"

    yield {
            "basename": "TASK: Compute linear model error from multiple models",
            "actions": [f"python {script}"],
            "file_dep": [script, script_linear_model, script_linear_model_error, exp_linear, exp_quadratic],
            "verbosity": 2, # show stdout
            }
