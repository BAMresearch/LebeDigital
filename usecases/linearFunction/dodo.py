from pathlib import Path
import yaml

def task_generate_virtual_samples():
    """
    Create virtual experiments to be stored in yaml files.
    """
    script = Path(__file__).parents[0] / "VirtualLinearModelExperiment.py"

    yield {
            "basename" : "TASK: generate virtual samples",
            "actions": ["python VirtualLinearModelExperiment.py"],
            "file_dep" : [script] ,
            "verbosity": 2, # show stdout
            }

def task_OptimizeLinearModel():
    """
    Optimize linear model based on experimental data stored in yaml files.
    """
    script = Path(__file__).parents[0] / "Optimize.py"
    script_linear_model = Path(__file__).parents[0] / "LinearModel.py"
    script_linear_model_error = Path(__file__).parents[0] / "LinearModelError.py"

    yield {
            "basename" : "TASK: Compute linear model error from multiple models",
            "actions": [f"python {script}"],
            "file_dep": [script, script_linear_model, script_linear_model_error],
            "verbosity": 2, # show stdout
            }
