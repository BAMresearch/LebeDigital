from pathlib import Path
import yaml

def task_generate_virtual_samples():
    """
    Create virtual experiments to be stored in yaml files.
    """
    script = Path(__file__).parents[0] / "virtual_linear_model_experiment.py"

    yield {
            "basename" : "TASK: generate virtual samples",
            "actions": ["python virtual_linear_model_experiment.py"],
            "file_dep" : [script] ,
            "verbosity": 2, # show stdout
            }

def task_OptimizeLinearModel():
    """
    Optimize linear model based on experimental data stored in yaml files.
    """
    script = Path(__file__).parents[0] / "optimize.py"
    script_linear_model = Path(__file__).parents[0] / "linear_model.py"
    script_linear_model_error = Path(__file__).parents[0] / "linear_model_error.py"

    yield {
            "basename" : "TASK: Compute linear model error from multiple models",
            "actions": [f"python {script}"],
            "file_dep": [script, script_linear_model, script_linear_model_error],
            "verbosity": 2, # show stdout
            }
