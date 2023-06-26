import os
from pathlib import Path

import numpy
import pytest

from lebedigital.calibration.utils import extract_third_load_cycle


def test_third_loading_cycle():
    # defining paths and directories
    data_dir = "calibration_data"
    data_path = Path(__file__).parent / data_dir
    input_file = "Wolf 8.2 Probe 1.csv"
    file_path = data_path / input_file

    extracted_data = extract_third_load_cycle(str(file_path), threshold=0.5)

    assert numpy.min(extracted_data["Force [kN]"]) == pytest.approx(-118.55971)
    assert numpy.max(extracted_data["Force [kN]"]) == pytest.approx(-6.1672144)
