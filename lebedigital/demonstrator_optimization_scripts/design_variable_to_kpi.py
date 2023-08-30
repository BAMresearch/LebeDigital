import json
import os, sys
from lebedigital.demonstrator_optimization_scripts.utils import load_json, update_json

def design_var_to_kpi(workflow_path:str,X: dict, seed: int) -> dict:
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
    # Pass the parameter to X to the input to forward. Meaning overwrrite the input.
    # The design variables, aggregate ratio and the slag ratio needs to be updated.

    #TODO: the below is hardcoded. fixit
    #design_var_paths = {'aggregates_volume_fraction': workflow_path +'/Inputs/aggregates_volume_fraction.json',
    #                    'sc_volume_fraction': workflow_path + '/Inputs/sc_fraction.json'}
    design_var_paths = {'height': workflow_path + '/Inputs/geometry.json',
                        'sc_mass_fraction': workflow_path + '/Inputs/sc_fraction.json'}

    for key, value in X.items():
        update_json(design_var_paths[key],key,value)

    # pass the seed to the scripts for the RVs (see eqn 29 SVO paper)
    # Updating the phi's which are input to the script.
    # phi_hydration_path = workflow_path + '/Inputs/phi_hydration.json'
    # phi_paste_path = workflow_path + '/Inputs/phi_paste.json'
    # update_json(phi_hydration_path, 'seed', seed)
    # update_json(phi_paste_path, 'seed', seed)

    seed_path = workflow_path + '/Inputs/seed_learnt_models.json'
    update_json(seed_path, 'seed', seed)

    # Run the workflow using snakemake
    # add the path to the workflow file and the path to the directory
    workflow_file_path = workflow_path + '/Snakefile'
    os.system(f'snakemake --cores 7 --snakefile {workflow_file_path} '
              f'--directory {workflow_path}  workflow_targets --use-conda')

    # Read in the KPIs in a dict
    Results_path = workflow_path + '/Results/'
    FEM_KPI = Results_path + 'kpi_from_fem.json'
    gwp_KPI = Results_path + 'gwp_beam.json'
    beam_design_KPI = Results_path + 'beam_design.json'
    y = {}
    for i, path in enumerate([FEM_KPI, gwp_KPI, beam_design_KPI]):
        tmp = load_json(path)
        y.update(tmp)

    # return the KPIs
    #TODO: this is specific to the constraints and objective choosen. careful
    kpi = {
        "gwp_beam": y["gwp_beam"]["value"],
       # "check_steel_area": y["check_steel_area"]["value"],
        "constraint_beam_design": y["constraint_beam_design"]["value"], 
        "constraint_temperature": y["constraint_temperature"]["value"],
        "constraint_time": y["constraint_time"]["value"]
    }
    return kpi

if __name__ == '__main__':
    path = '../../usecases/optimization_paper/1'
    design_var = {'height': 260,
                        'sc_mass_fraction': 0.35}
    seed = 66
    design_var_to_kpi(workflow_path=path,X=design_var,seed=seed)