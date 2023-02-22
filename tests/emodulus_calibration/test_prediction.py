from lebedigital.calibration.posterior_predictive import perform_prediction_three_point_bending
import pytest
import numpy as np
import os
from pathlib import Path


def test_prediction():
    print("Prediction test started")

    # check if KG is present
    kg_path = "prediction_data/prediction_kg.owl"
    if os.path.exists(os.path.join(os.path.dirname(Path(__file__)),kg_path)):
        path_calibrated_data = os.path.join(os.path.dirname(Path(__file__)),kg_path)
    else:
        print("The knowledge graph file is not present")

    pos_pred = perform_prediction_three_point_bending(knowledge_graph_file=path_calibrated_data, mode="cheap")
    print(f"The posterior predictive mean is {np.mean(pos_pred)} N/mm2")
    assert np.mean(pos_pred) == pytest.approx(125.00497593013911)  # this needs to specified by subject matter experts

    print("Prediction test completed")
