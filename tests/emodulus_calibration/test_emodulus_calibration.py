from lebedigital.calibration.calibrationWorkflow import estimate_youngs_modulus
from lebedigital.calibration.utils import read_exp_data_E_mod
import pytest
import numpy as np
from pathlib import Path

def test_emodulus_calibration():
    # defining paths and directories
    data_dir = 'calibration_data'
    data_path = Path(__file__).parent / data_dir
    input_file = "Wolf 8.2 Probe 1.csv"

    # defining calibration input, setting the values for the priors
    E_loc = 30  # KN/mm2 (mean)
    E_scale = 10  # KN/mm2 (std)

    # query knowledge graph
    # Note remove the mode 'dheap' for local testing so that 'full' sampling can be tested
    # passing the length and diameter for the specific experiment
    output = read_exp_data_E_mod(path=data_path, exp_name=input_file,length=300.2,diameter=98.6)
    E_samples = estimate_youngs_modulus(experimental_data=output,
                                        calibration_metadata={"E_loc": E_loc, "E_scale": E_scale},
                                        calibrated_data_path=data_path, mode='cheap')

    # checking if the KG file for the calibration is created
    assert (data_path / f'calibrationWorkflow{input_file}.owl').is_file()

    # checking for the mean of the calibrated E modulus
    # Note : Somehow the seed in EMCEE deosnt seed to work.
    assert np.mean(E_samples) == pytest.approx(31, rel=0.2)


