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
    E_samples, post_pred = esimate_Youngs_modulus(
        experimental_data=output, calibration_metadata={"E_loc": 100, "E_scale": 100},
        calibrated_data_path=os.path.dirname(Path(__file__)), mode='cheap'
    )
    # checking if the KG file for the calibration is created
    assert os.path.exists(os.path.join(os.path.dirname(Path(__file__)), 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl')) == True

    # checking for the calibrated E modulus and the posterior predictive
    assert np.mean(post_pred) == pytest.approx(95, rel=50)

    assert np.mean(E_samples) == pytest.approx(95E03, rel=50e03)

    print("Calibration test completed")

assert os.path.exists(os.path.join(os.path.dirname(Path(__file__)), 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl')) == True
# if __name__ == "__main__":
#     output = query_KG(
#         path=os.path.dirname(Path(__file__)), exp_name="Wolf 8.2 Probe 1.csv"
#     )
#     E_samples, post_pred = esimate_Youngs_modulus(
#         experimental_data=output, calibration_metadata={"E_loc": 100, "E_scale": 100},
#         calibrated_data_path=os.path.dirname(Path(__file__)),mode='cheap'
#     )
#

