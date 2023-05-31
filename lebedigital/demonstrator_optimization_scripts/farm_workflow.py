import os
import shutil
import numpy as np
import sys

def farm_workflow(path:str,seed:list):
    """
    Created multiple copies of snakemake workflow folder to be used later to parallelize
    Parameters
    ----------
    path: should point to the optimization_paper path
    seed

    Returns
    -------

    """

    for i,v in enumerate(seed):
        new_dir_path = os.path.join(path,str(i+1))
        src_path = path + '/optimization_workflow'
        # copy to the newly created folder
        if not os.path.exists(new_dir_path):
            shutil.copytree(src=src_path,dst=new_dir_path)
        else:
            shutil.rmtree(new_dir_path)
            shutil.copytree(src=src_path,dst=new_dir_path)

if __name__ == '__main__':
    path = '../../usecases/optimization_paper'
    #seed = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    seed = np.load('../../usecases/demonstrator/Calibration/seed_tmp.npy').astype(int).tolist()
    farm_workflow(path,seed)
