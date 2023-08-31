import json
from pathlib import Path

import yaml


def py_macros(file_list):
    for file_name in file_list:
        file_name = str(file_name)  # convert pathlib object to useful string

        # read dictionaries
        with open(file_name + ".yaml") as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)

        # write tex macros
        with open(file_name + ".tex", "w+") as f:
            for dictionary in data:
                for key in data[dictionary]:
                    f.write(f"\\newcommand{{\\{key}}}{{{data[dictionary][key]}}}\n")


def input_optimization_macros(path_to_workflow, path_to_tex_file):
    # get a list of files in path to workflow directory
    import os

    with open(path_to_tex_file, "w+") as f:
        f.write("% This file is generated automatically by py_macros.py\n")

    input_path = Path(path_to_workflow) / "Inputs"
    input_files = os.listdir(input_path)

    replace_dict = {
        "28": "twentyeight",
        "_": "",
        "0": "zero",
        "1": "one",
        "2": "two",
        "3": "three",
        "4": "four",
        "5": "five",
        "6": "six",
        "7": "seven",
        "8": "eight",
        "9": "nine",
    }

    for input_file in input_files:
        # open file, update dictionary, write file
        with open(input_path / input_file) as f:
            file_data = yaml.load(f, Loader=yaml.SafeLoader)

        # create tex macros
        with open(path_to_tex_file, "a") as f:
            for variable in file_data:
                variable_name = variable
                for word, replacement in replace_dict.items():
                    variable_name = variable_name.replace(word, replacement)

                variable_name = "input" + variable_name
                f.write(f"\\newcommand{{\\{variable_name}}}{{{file_data[variable]['value']}}}\n")
                f.write(f"\\newcommand{{\\{variable_name + 'unit'}}}{{{file_data[variable]['unit']}}}\n")

    result_path = Path(path_to_workflow) / "Results"

    # check if results directory exists
    if not result_path.exists():
        # create directory
        os.mkdir(result_path)

    # get a list with all json files in result_path directory
    result_files = [f for f in os.listdir(result_path) if f.endswith(".json")]

    for result_file in result_files:
        # open file, update dictionary, write file
        with open(result_path / result_file) as f:
            file_data = yaml.load(f, Loader=yaml.SafeLoader)
        #
        # # create tex macros
        with open(path_to_tex_file, "a") as f:
            for variable in file_data:
                variable_name = variable
                for word, replacement in replace_dict.items():
                    variable_name = variable_name.replace(word, replacement)

                variable_name = "results" + variable_name
                f.write(f"\\newcommand{{\\{variable_name}}}{{{file_data[variable]['value']}}}\n")
                f.write(f"\\newcommand{{\\{variable_name + 'unit'}}}{{{file_data[variable]['unit']}}}\n")
