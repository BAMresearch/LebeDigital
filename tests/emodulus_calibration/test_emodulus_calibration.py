from lebedigital.calibration.calibrationWorkflow import esimate_Youngs_modulus
from lebedigital.calibration.utils import query_KG
import pytest
import numpy as np
from pathlib import Path

def test_emodulus_calibration():
    # defining paths and directories
    data_dir = 'calibration_data'
    data_path = Path(__file__).parent / data_dir
    input_file = "Wolf 8.2 Probe 1.csv"

    # defining calibration input
    E_loc = 30  # TODO (Atul) what is this??? this needs units!
    E_scale = 10  # TODO (Atul) what is this???

    # query knowledge graph
    output = query_KG(path=data_path, exp_name=input_file)
    E_samples = esimate_Youngs_modulus(experimental_data=output,
                                       calibration_metadata={"E_loc": E_loc, "E_scale": E_scale},
                                       calibrated_data_path=data_path, mode='cheap')

    # checking if the KG file for the calibration is created
    assert (data_path / f'calibrationWorkflow{input_file}.owl').is_file()

    # checking for the mean of the calibrated E modulus
    # TODO (Atul), there should be a way to make sure the result is the same each time? Defining a random seed?
    #      rel = 0.4 is too large...
    assert np.mean(E_samples) == pytest.approx(E_loc, rel=0.4)


