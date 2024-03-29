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
    mix_gwp = load_json(results_path / "gwp_mix.json")
    kpis.update(mix_gwp)
    gwp_steel = load_json(results_path / "steel_gwp_per_volume.json")
    kpis.update(gwp_steel)

    # return the KPIs
    return kpis


if __name__ == "__main__":
    path_to_workflow = Path("../optimization_workflow")
    input_path = path_to_workflow / "Inputs"
    output_path = path_to_workflow / "Results"

    # input lists
    height_list = [1200]
    slag_ratio_list = [0.5]

    df = pd.DataFrame()

    for i, height in enumerate(height_list):
        for j, slag_ratio in enumerate(slag_ratio_list):
            # remove all files from the directory output_path
            for file in os.listdir(output_path):
                os.remove(output_path / file)

            total = len(height_list) * len(slag_ratio_list)
            current = i * len(slag_ratio_list) + j + 1
            print("___________________________________________________________")
            print(f" {current}/{total}     RUN WORKFLOW WITH {height} {slag_ratio}")
            print("___________________________________________________________")
            inputs = {"height": height, "slag_ratio": slag_ratio}
            results = get_kpis(inputs, path_to_workflow)

            new_row = {
                "height": inputs["height"],
                "slag_ratio": inputs["slag_ratio"],
                "gwp": results["gwp_beam"]["value"],
                "constraint_beam_design": results["constraint_beam_design"]["value"],
                "constraint_temperature": results["constraint_temperature"]["value"],
                "constraint_time": results["constraint_time"]["value"],
            }

            # build new pandas dataframe from dictionary
            new_df = pd.DataFrame(new_row, index=[0])
            # add new row to existing dataframe
            df = pd.concat([df, new_df], ignore_index=True)

    # df.to_csv(f"kpis_{inputs['agg_ratio']}_{inputs['slag_ratio']}.csv",index=False)
    df.to_csv(f"kpis.csv", index=False)

print("Done")
