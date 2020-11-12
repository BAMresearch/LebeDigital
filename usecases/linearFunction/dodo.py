from pathlib import Path
import yaml

def task_computepyMC3():
    """
    Create one task for each yaml file.
    """
    parameter_files = (Path(__file__).parents[0]).glob('*.yaml')
    script = Path(__file__).parents[0] / "ComputepyMC3.py"
    script_general = Path(__file__).parents[0] / "Correlation.py"

    for parameter_file in parameter_files:
        pickle_file = parameter_file.with_suffix(".pkl")
        action = f"python {script} {parameter_file}"
        yield {
                "basename" : "TASK: compute pyMC3",
                "name" : parameter_file.name,
                "actions": [action],
                "file_dep" : [script, script_general, parameter_file] ,
                "targets" : [pickle_file],
                "verbosity": 2, # show stdout
                }

def task_PostProcess():
    """
    Create one task for each pickle (trace) file.
    """
    parameter_files = (Path(__file__).parents[0]).glob('*.yaml')
    script = Path(__file__).parents[0] / "PostProcess.py"
    script_general = Path(__file__).parents[0] / "Correlation.py"

    for parameter_file in parameter_files:
        plot_file = parameter_file.with_suffix(".png")
        pickle_file = parameter_file.with_suffix(".pkl")

        action = f"python {script} {parameter_file}"
        yield {
                "basename" : "TASK: PostProcess",
                "name" : parameter_file.name,
                "actions": [action],
                "file_dep" : [script, script_general, pickle_file] ,
                "targets" : [plot_file],
                "verbosity": 2, # show stdout
                }

