from pathlib import Path
import json
import pandas as pd
import numpy as np

import os

def update_json(file_path: Path, key: str, value):
    # Read the JSON file
    with open(file_path, 'r') as f:
        data = json.load(f)
    # TODO:will work only when 'value' key is present
    # Update the value of the specified key
    data[key]['value'] = value

    # Write the updated data back to the JSON file
    with open(file_path, 'w') as f:
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
    variables. The Random variables (b) x->b->KPIs are also sampled for a given value of seed.
    Args:
     X: (dict) with keys 'agg_ratio' (volume ratio of the aggregates) and 'slag_ratio'
     seed: the seed parameter. This ensures that the sampled Random variable here is the same as the one passed in the
     forward call

    Returns:
        y : dict with all the KPIs

    """

    input_path = path / 'Inputs'
    # Pass the parameter to X to the input to forward. Meaning overwrrite the input.
    # The design variables, aggregate ratio and the slag ratio needs to be updated.
    update_json(input_path / 'aggregates_volume_fraction.json', 'aggregates_volume_fraction', input['agg_ratio'])
    update_json(input_path / 'sc_fraction.json', 'sc_volume_fraction', input['slag_ratio'])

    # # pass the seed to the scripts for the RVs (see eqn 29 SVO paper)
    # # Updating the phi's which are input to the script.
    # update_json(phi_hydration_path, 'seed', seed)
    # update_json(phi_paste_path, 'seed', seed)

    # Run the workflow using snakemake
    # add the path to the workflow file and the path to the directory
    # workflow_file_path = Optimization_workflow_path + '/Snakefile'
    # directory_path = Optimization_workflow_path

    # run workflow
    os.system(f'snakemake --cores 4 --snakefile {path / "Snakefile"} '
              f'--directory {path}')


    # get kpis
    kpis = {}
    results_path = path / 'Results'
    kpi_from_fem = load_json(results_path / 'kpi_from_fem.json')
    kpis.update(kpi_from_fem)
    gwp_beam = load_json(results_path / 'gwp_beam.json')
    kpis.update(gwp_beam)
    beam_design = load_json(results_path / 'beam_design.json')
    kpis.update(beam_design)

    # return the KPIs
    return kpis


if __name__ == "__main__":

    path_to_workflow = Path('../optimization_workflow')
    input_path = path_to_workflow / 'Inputs'

    # input lists
    aggregate_ratio_list = [0.1, 0.5, 0.8]
    slag_ratio_list = [0.1, 0.5, 0.8]

    df = pd.DataFrame()

    for i, agg_ratio in enumerate(aggregate_ratio_list):
        for j, slag_ratio in enumerate(slag_ratio_list):
            total = len(aggregate_ratio_list) * len(slag_ratio_list)
            current = i*len(slag_ratio_list) + j +1
            print('___________________________________________________________')
            print(f' {current}/{total}     RUN WORKFLOW WITH {agg_ratio} {slag_ratio}')
            print('___________________________________________________________')
            inputs = {'agg_ratio': agg_ratio, 'slag_ratio': slag_ratio}
            results = get_kpis(inputs, path_to_workflow)

            df = df.append({'agg_ratio': inputs['agg_ratio'],
                                 'slag_ratio': inputs['slag_ratio'],
                                 'gwp': results['gwp_mix']['value'],
                                 'check_beam_design': results['check_steel_area']['value'],
                                 'max_temp': results['max_reached_temperature']['value'],
                                 'time_of_demoulding': results['time_of_demolding']['value']
                                 }, ignore_index=True)

    #df.to_csv(f"kpis_{inputs['agg_ratio']}_{inputs['slag_ratio']}.csv",index=False)
    df.to_csv(f"kpis.csv",index=False)

print('Done')