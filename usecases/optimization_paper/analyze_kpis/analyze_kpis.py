import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def update_json(file_path: Path, key: str, value):
    # Read the JSON file
    with open(file_path, "r") as f:
        data = json.load(f)
    # TODO:will work only when 'value' key is present
    # Update the value of the specified key
    data[key]["value"] = value

    # Write the updated data back to the JSON file
    with open(file_path, "w") as f:
        json.dump(data, f, indent=4, sort_keys=True)


def load_json(path: Path) -> dict:
    data = None
    if path.is_file():
        with open(path) as f:
            data = json.load(f)
    return data


def get_kpis(input: dict, path: Path) -> dict:
    """
    Runs the snakemake workflow and the returns the KPIs for objective and constraints for a given value of the design
    variables.
    Args:
        input: dict with the design variables and the aggregate ratio and slag ratio
        path: Path to the workflow directory
    Returns:
        kpis : dict with all the KPIs

    """
    input_path = path / "Inputs"
    # Pass the parameter to X to the input to forward. Meaning overwrite the input.
    # The design variables, aggregate ratio and the slag ratio needs to be updated.
    update_json(input_path / "geometry.json", "height", input["height"])
    update_json(input_path / "sc_fraction.json", "sc_mass_fraction", input["slag_ratio"])

    # # pass the seed to the scripts for the RVs (see eqn 29 SVO paper)
    # # Updating the phi's which are input to the script.
    # update_json(phi_hydration_path, 'seed', seed)
    # update_json(phi_paste_path, 'seed', seed)

    # Run the workflow using snakemake
    # add the path to the workflow file and the path to the directory
    # workflow_file_path = Optimization_workflow_path + '/Snakefile'
    # directory_path = Optimization_workflow_path

    # run workflow
    os.system(f'snakemake --cores 4 --snakefile {path / "Snakefile"} ' f"--directory {path}")

    # get kpis
    kpis = {}
    results_path = path / "Results"
    kpi_from_fem = load_json(results_path / "kpi_from_fem.json")
    kpis.update(kpi_from_fem)
    gwp_beam = load_json(results_path / "gwp_beam.json")
    kpis.update(gwp_beam)
    beam_design = load_json(results_path / "beam_design.json")
    kpis.update(beam_design)

    # return the KPIs
    return kpis


if __name__ == "__main__":
    path_to_workflow = Path("../optimization_workflow")
    input_path = path_to_workflow / "Inputs"

    # input lists
    height_list = [500.0, 700.0, 900.0]
    slag_ratio_list = [0.1, 0.5, 0.8]

    df = pd.DataFrame()

    data = {
        "height": [],
        "slag_ratio": [],
        "gwp": [],
        "check_beam_design": [],
        "max_temp": [],
        "time_of_demoulding": [],
    }

    for i, height in enumerate(height_list):
        for j, slag_ratio in enumerate(slag_ratio_list):
            total = len(height_list) * len(slag_ratio_list)
            current = i * len(slag_ratio_list) + j + 1
            print("___________________________________________________________")
            print(f" {current}/{total}     RUN WORKFLOW WITH {height} {slag_ratio}")
            print("___________________________________________________________")
            inputs = {"height": height, "slag_ratio": slag_ratio}
            results = get_kpis(inputs, path_to_workflow)

            data["height"].append(inputs["height"])
            data["slag_ratio"].append(inputs["slag_ratio"])
            data["gwp"].append(results["gwp_beam"]["value"])
            data["check_beam_design"].append(results["constraint_beam_design"]["value"])
            data["max_temp"].append(results["max_reached_temperature"]["value"])
            data["time_of_demoulding"].append(results["time_of_demolding"]["value"])

    df = pd.DataFrame(data)
    df.to_csv(f"kpis.csv", index=False)

print("Done")
