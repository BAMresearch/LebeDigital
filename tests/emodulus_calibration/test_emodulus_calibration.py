from lebedigital.calibration.calibrationWorkflow import esimate_Youngs_modulus
from lebedigital.calibration.utils import query_KG
from probeye.ontology.knowledge_graph_import import import_parameter_samples

import pytest
import numpy as np
import os
from pathlib import Path

def test_emodulus_calibration():

    # query knowledge graph. Note KG team needs to provide this.
    print("Calibration test started")
    output = query_KG(
        path=os.path.dirname(Path(__file__)), exp_name="Wolf 8.2 Probe 1.csv"
    )
    E_samples= esimate_Youngs_modulus(
        experimental_data=output, calibration_metadata={"E_loc": 30, "E_scale": 10},
        calibrated_data_path=os.path.dirname(Path(__file__)), mode='cheap'
    )
    # checking if the KG file for the calibration is created
    assert os.path.exists(os.path.join(os.path.dirname(Path(__file__)), 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl')) == True

    # checking for the calibrated E modulus
    print(f"The posterior E mean is {np.mean(E_samples)} KN/mm2")
    assert np.mean(E_samples) == pytest.approx(30, rel=0.3) # this needs to specified by subject matter experts

    print("Calibration test completed")


