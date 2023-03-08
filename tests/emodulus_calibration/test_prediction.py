from lebedigital.calibration.posterior_predictive import perform_prediction_three_point_bending
import pytest
import numpy as np
import os
from pathlib import Path


def test_prediction():
    # defining paths and directories
    # TODO: this is not good!!! one test depending on another
    #       however, as we are planning on removing the ontologie part anyways, this will soon change...
    data_dir = 'calibration_data'
    data_path = Path(__file__).parent / data_dir
    input_file = 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl'
    kg_path = data_path / input_file

    # check if KG is present
    assert kg_path.is_file()

    pos_pred = perform_prediction_three_point_bending(knowledge_graph_file=str(kg_path), mode="cheap")

    # TODO carry units till the end (if possible...)!!!
    #      this result depends on the result of the calibration, which is random...
    #      this needs to change once we fix the interfaces and datatypes..
    assert np.mean(pos_pred) == pytest.approx(125, rel=0.3)  # rel=0.5, only for debugging
