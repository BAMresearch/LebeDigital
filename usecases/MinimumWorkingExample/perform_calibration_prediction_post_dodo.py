import os
from pathlib import Path
import numpy as np
import pandas as pd

import yaml

from lebedigital.calibration.calibrationWorkflow import estimate_youngs_modulus
from lebedigital.calibration.utils import read_exp_data_E_mod
from lebedigital.calibration.posterior_predictive_three_point_bending import perform_prediction
from lebedigital.calibration.posterior_predictive_three_point_bending import wrapper_three_point_bending

# when "cheap option" or "single" is run, only this souce of raw data is processed
single_example_name = 'Wolf 8.2 Probe 1'

# corresponding mix sor the "single" example
single_mix_name = '2014_12_10 Wolf.xls'

ParentDir = os.path.dirname(Path(__file__))
emodul_output_directory = Path(ParentDir, 'emodul')

metadata_emodulus_directory = Path(
    emodul_output_directory, 'metadata_yaml_files')
processed_data_emodulus_directory = Path(
    emodul_output_directory, 'processed_data')  # folder with csv data files

calibrated_data_directory = Path(emodul_output_directory, 'calibrated_data')
predicted_data_directory = Path(emodul_output_directory,'predicted_data')

config = {'mode': 'full'}
def task_perform_calibration():
    """Loop over the experiments and store the calibrated E in csv file. Each iteration generates four files
    viz. calibrated data.csv, displacement data used, load data used and the knowledge graph file containing the calibration details
    (these are not specific file names)"""
    # create folder, if it is not there
    Path(calibrated_data_directory).mkdir(parents=True, exist_ok=True)

    # defining calibration input, setting the values for the priors
    E_loc = 30  # KN/mm2 (mean)
    E_scale = 10  # KN/mm2 (std)

    # setting for fast test, defining the list
    if config['mode'] == 'cheap':
        list_exp_name = [single_example_name]
    else:  # go through all files
        list_exp_name = os.listdir(
            processed_data_emodulus_directory)
        list_exp_name = [os.path.splitext(f)[0] for f in list_exp_name] # split extension
    for f in list_exp_name:
        exp_metadata_path = os.path.join(metadata_emodulus_directory,f+'.yaml')
        with open(exp_metadata_path) as file:
            data = yaml.safe_load(file)
        diameter = float(data['diameter'])
        length = float(data['length'])
        output = read_exp_data_E_mod(path=processed_data_emodulus_directory, exp_name=f+'.csv',
                                     length=length,diameter=diameter)
        E_samples = estimate_youngs_modulus(experimental_data=output,
                                            calibration_metadata={"E_loc": E_loc, "E_scale": E_scale},
                                            calibrated_data_path=calibrated_data_directory, mode=config['mode'])

        # store samples as csv
        np.savetxt(os.path.join(calibrated_data_directory,f+'_calibrated_samples.csv'),E_samples,delimiter=',')


def task_perform_prediction():
    # create folder, if it is not there
    Path(predicted_data_directory).mkdir(parents=True, exist_ok=True)

    # choose calibrated data to be used for prediction
    calibrated_data = single_example_name
    df = pd.read_csv(os.path.join(calibrated_data_directory,single_example_name+'_calibrated_samples.csv'), header=None)
    samples = df[0].tolist()

    # perform prediction
    pos_pred = perform_prediction(forward_solver=wrapper_three_point_bending,parameter=samples)

    # store prediction in csv
    np.savetxt(os.path.join(predicted_data_directory,'prediction_with_'+calibrated_data+'.csv'),pos_pred,delimiter=',')

if __name__ == '__main__':
    task_perform_calibration()
    task_perform_prediction()

