from lebedigital.calibration.posterior_predictive_three_point_bending import perform_prediction_three_point_bending
from lebedigital.calibration.posterior_predictive_three_point_bending import wrapper_three_point_bending
from probeye.ontology.knowledge_graph_import import import_parameter_samples
import pytest
import numpy as np
import os
from pathlib import Path


def test_prediction():
    # Reading in the calibrated data. TO be pas here by some KG
    E_samples = [30.1,30.2,30.53,30.8,29.6] #kN/mm2
    KG = False
    if KG:
        data_dir = 'calibration_data'
        data_path = Path(__file__).parent / data_dir
        input_file = 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl'
        kg_path = data_path / input_file
        assert kg_path.is_file()
        E_samples = import_parameter_samples(kg_path)

    # performing posterior predictive
    pos_pred = perform_prediction_three_point_bending(forward_solver=wrapper_three_point_bending,
                                                      parameter=E_samples)
    assert np.mean(pos_pred) == pytest.approx(120, rel=0.1)  # rel=0.5, only for debugging
