import os
from pathlib import Path

import numpy as np
import pytest

from lebedigital.calibration.calibrationWorkflow import (
    _check_E_mod_calibration_metadata, _check_E_mod_experimental_data,
    estimate_youngs_modulus)
from lebedigital.calibration.utils import read_exp_data_E_mod

# for reproducible results
seed =1
np.random.seed(seed)

def test_emodulus_calibration():
    # defining paths and directories
    data_dir = "calibration_data"
    data_path = Path(__file__).parent / data_dir
    input_file = "Wolf 8.2 Probe 1.csv"

    # defining calibration input, setting the values for the priors
    E_loc = 30  # KN/mm2 (mean)
    E_scale = 10  # KN/mm2 (std)

    # query knowledge graph
    # Note remove the mode 'dheap' for local testing so that 'full' sampling can be tested
    # passing the length and diameter for the specific experiment
    output = read_exp_data_E_mod(path=data_path, exp_name=input_file, length=300.2, diameter=98.6)
    E_samples = estimate_youngs_modulus(
        experimental_data=output,
        calibration_metadata={"E_loc": E_loc, "E_scale": E_scale},
        calibrated_data_path=data_path,
        mode="test",
    )

    # checking if the KG file for the calibration is created
    assert (data_path / f"calibrationWorkflow{os.path.splitext(input_file)[0]}").is_file()

    # checking for the mean of the calibrated E modulus
    assert np.mean(E_samples) == pytest.approx(31.334306458960317)


def test_check_E_mod_calibration_metadata():
    calibration_metadata = {"E_loc": 30, "E_scale": 10}
    assert _check_E_mod_calibration_metadata(calibration_metadata) == True

    calibration_metadata = {"E_loc": 30}
    assert _check_E_mod_calibration_metadata(calibration_metadata) == False


def test_check_E_mod_experimental_data():
    experimental_data = {
        "exp_name": "Wolf 8.2 Probe 1.csv",
        "force": [1, 2, 3],
        "displacement": [1, 2, 3],
        "height": 300.2,
        "diameter": 98.6,
    }
    assert _check_E_mod_experimental_data(experimental_data) == True

    experimental_data = {
        "exp_name": "Wolf 8.2 Probe 1.csv",
        "force": [1, 2, 3],
        "displacement": [1, 2, 3],
        "height": 300.2,
    }
    assert _check_E_mod_experimental_data(experimental_data) == False
