from lebedigital.calibration.calibrationWorkflow import perform_calibration

import pytest
import numpy as np
import os
from pathlib import Path

def test_emodulus_calibration():
    ParentDir = os.path.dirname(os.path.dirname(os.path.dirname(Path(__file__))))
    path_KG = "../../usecases/MinimumWorkingExample/emodul/knowledge_graphs"  # This is not used so doesnt matter.
    path_calibrated = ParentDir + "/usecases/MinimumWorkingExample/emodul/calibrated_data"
    # TODO: add a seed in the calibration
    E_calibrated = perform_calibration(path_to_KG=path_KG,calibrated_data_path=path_calibrated,mode='cheap',KnowledgeGraph=False)
    assert np.mean(E_calibrated) == pytest.approx(95E03,rel = 50e03) # Just adding test for calibration, not for prediction
