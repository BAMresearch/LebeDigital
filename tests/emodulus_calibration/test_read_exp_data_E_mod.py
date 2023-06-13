from lebedigital.calibration.utils import read_exp_data_E_mod
import pytest
import numpy as np
from pathlib import Path

def test_read_exp_data_E_mod():
    # defining paths and directories
    data_dir = 'calibration_data'
    data_path = Path(__file__).parent / data_dir
    input_file = "Wolf 8.2 Probe 1.csv"

    # the length, diameter to be provided by the KG module
    output = read_exp_data_E_mod(path=data_path, exp_name=input_file,length=300.2,diameter=98.6)

    assert np.min(output["force"]) == pytest.approx(-118.55971)
    assert np.max(output["force"]) == pytest.approx(-6.1672144)
    assert output["diameter"] == pytest.approx(98.6)