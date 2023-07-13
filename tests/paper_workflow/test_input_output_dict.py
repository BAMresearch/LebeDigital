import json
import os
from pathlib import Path


def test_input_output_dict():
    """
    this test makes sure that the keys in the input and output dictionaries are consistent
    this is especially important for the workflow and for the tex macros
    the test only makes sense once the snakemake workflow has been run once
    """

    # assuming the tests are started in the tests folder
    # path to current working directory
    path_to_cwd = Path(os.getcwd())
    path_to_base = path_to_cwd

    while path_to_base.name != "LebeDigital":
        path_to_base = path_to_base.parent

    path_to_workflow = path_to_base / Path("usecases/optimization_paper/optimization_workflow/")
    path_to_input = path_to_workflow / "Inputs"
    path_to_output = path_to_workflow / "Results"

    file_list = []

    if path_to_input.is_dir():
        file_list += [path_to_input / f for f in os.listdir(path_to_input) if f.endswith(".json")]
    if path_to_output.is_dir():
        file_list += [path_to_output / f for f in os.listdir(path_to_output) if f.endswith(".json")]

    dict = {}
    for path in file_list:
        with open(path) as f:
            content = json.loads(f.read())
            for key in content.keys():
                assert key not in dict.keys()

            dict.update(content)
