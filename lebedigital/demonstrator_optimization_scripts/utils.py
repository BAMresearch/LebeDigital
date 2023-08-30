import os, sys
import json

import numpy as np

from lebedigital.demonstrator_optimization_scripts.farm_workflow import farm_workflow

def python_fn_run_jobs(path_to_scripts_folder:str,no_samples:int):
    """
    Create copies of the workflow for each sample, and run slurm job array.
    Parameters
    ----------
    path_to_scripts_folder
    no_samples

    Returns
    -------

    """
    #current_dir = os.path.dirname(os.path.abspath(__file__))
    #path_farm = path_to_scripts_folder + "/farm_workflow.py"
    #path_jobs = path_to_scripts_folder + 'run_jobs.sh'

    print("!!! Creating folders for copies of the workflow !!!")
    #os.system(f'python {path_farm}')
    seed = np.load(path_to_scripts_folder+'../../usecases/demonstrator/Calibration/seed_tmp.npy').astype(int).tolist()
    farm_workflow(path=path_to_scripts_folder+'../../usecases/optimization_paper',
                  seed=seed)

    print("!!! folder creating DONE !!!")
    print("!!! Run multiple jobs in cluster !!!")
    # change directory to the file in which this fn is.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    original_dir = os.getcwd()
    os.chdir(script_dir)
    os.system(f'sbatch --wait --array=1-{no_samples} run_jobs.sh')
    if not os.path.exists(f'../../usecases/optimization_paper/{no_samples}/kpi.json'):
        raise FileNotFoundError
    print('All jobs finished')
    # restore to the working directory
    os.chdir(original_dir)

# %%
def load_json(path: str) -> dict:
    if path[-5:] == '.json':
        with open(path) as f:
            data = json.load(f)
    return data


# %%
def update_json(file_path: str, key: str, value):
    # Read the JSON file
    with open(file_path, 'r') as f:
        data = json.load(f)
    # TODO:will work only when 'value' key is present
    # Update the value of the specified key
    data[key]['value'] = value

    # Write the updated data back to the JSON file
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4, sort_keys=True)

def read_kpis(kpi_path:str):
    """
    Read in the kpis from the
    Parameters
    ----------
    kpi_path

    Returns
    -------

    """

    if not os.path.exists(kpi_path):
        print(f"Error: File {kpi_path} does not exist.")
    data = load_json(kpi_path)
    #TODO: the below is specific to the problem
    # print("!!! Attention the KPIs are specific and can change. Careful.")
    obj = data["gwp_beam"]
    C_1 = data["constraint_beam_design"]
    C_2 = data["constraint_temperature"]
    C_3 = data["constraint_time"]

    return obj, C_1, C_2, C_3

if __name__=='__main__':
    #o,c1,c2,c3 =read_kpis(kpi_path='../../usecases/optimization_paper/1/kpi.json')
    #print(o,c1,c2,c3)

    # test function
    python_fn_run_jobs('./',5)


