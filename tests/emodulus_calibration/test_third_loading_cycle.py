from lebedigital.calibration.utils import extract_third_load_cycle

import numpy
import pytest

def test_third_loading_cycle():
    path_data = "/home/atul_0711/Documents/PhD_Tasks/LeBeDigital/Codes/ModelCalibration/ModelCalibration/tests/emodulus_calibration/Wolf 8.2 Probe 1.csv"
    extracted_data = extract_third_load_cycle(path_data, threshold=0.5, vizualize=True)

    assert numpy.min(extracted_data["Force [kN]"]) == pytest.approx(-118.55971)
    assert numpy.max(extracted_data["Force [kN]"]) == pytest.approx(-6.1672144)

