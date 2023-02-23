from lebedigital.calibration.posterior_predictive import perform_prediction_three_point_bending
import pytest
import numpy as np
import os
from pathlib import Path


def test_prediction():
    # defining paths and directories
    data_dir = 'prediction_data'
    data_path = Path(__file__).parent / data_dir
    input_file = 'calibrationWorkflowWolf 8.2 Probe 1.csv.owl'
    kg_path = data_path / input_file

    # check if KG is present
    if not kg_path.is_file():
        raise Exception('The knowledge graph file is not present')

    pos_pred = perform_prediction_three_point_bending(knowledge_graph_file=str(kg_path), mode="cheap")

    # TODO carry units till the end (if possible...)!!!
    assert np.mean(pos_pred) == pytest.approx(119.66482703410645)
