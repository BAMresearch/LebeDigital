import json
import sys
import numpy as np
from lebedigital.demonstrator_optimization_scripts.design_variable_to_kpi import design_var_to_kpi

def parallel_workflow(array_id,design_var:np.ndarray, seed:list):
    """
    For the given array_id of job, points the created folder, runs the workflow and
    saves a kpi dictionary.
    Parameters
    ----------
    array_id
    design_var: rows are the samples and columns the design variable number.
    seed

    Returns
    -------

    """
    path = '../../usecases/optimization_paper/' + str(array_id)
    # TODO: pass the below aslo, not hardcode

    idx = array_id -1
    # design_var = {'aggregates_volume_fraction':design_var[idx,0], #0.4
    #               'sc_volume_fraction': design_var[idx,1]} #0.35
    design_var = {'height': design_var[idx, 0],  # 0.4
                  'sc_volume_fraction': design_var[idx, 1]}  # 0.35
    kpi = design_var_to_kpi(workflow_path=path, X=design_var, seed=seed[idx])
    kpi_path = path + '/kpi.json'
    with open(kpi_path, 'w') as f:
        json.dump(kpi, f, indent=4, sort_keys=True)

# to pass the job array number here.
#TODO: read from the seed file
#seeds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

#seed = np.random.randint(666,size=5)
#np.save('../../usecases/demonstrator/Calibration/seed_tmp.npy',seed)
#des_var = np.random.uniform(size=(5,2))
#np.save('../../usecases/demonstrator/Calibration/design_var_tmp.npy',des_var)

seeds = np.load('../../usecases/demonstrator/Calibration/seed_tmp.npy').astype(int).tolist()
design_var = np.load('../../usecases/demonstrator/Calibration/design_var_tmp.npy')
parallel_workflow(int(sys.argv[1]),design_var=design_var,seed=seeds)
#parallel_workflow(4,design_var=design_var,seed=seeds)