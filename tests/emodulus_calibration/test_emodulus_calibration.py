from lebedigital.calibration.calibrationWorkflow import perform_calibration

import pytest
import numpy as np
import os
from pathlib import Path

def test_emodulus_calibration():
    ParentDir = os.path.dirname(os.path.dirname(os.path.dirname(Path(__file__))))
    path_KG = ParentDir + "/tests/emodulus_calibration/Wolf 8.2 Probe 1.csv"  # This is not used so doesnt matter.
    path_calibrated = ParentDir + "/tests/emodulus_calibration"
    print("Calibration test started")
    # TODO: add a seed in the calibration
    E_calibrated = perform_calibration(path_to_KG=path_KG,calibrated_data_path=path_calibrated,mode='cheap',KnowledgeGraph=False)
    assert np.mean(E_calibrated) == pytest.approx(95E03,rel = 50e03) # Just adding test for calibration, not for prediction
    print("Calibration test completed")