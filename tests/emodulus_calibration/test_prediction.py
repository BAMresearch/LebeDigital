from lebedigital.calibration.posterior_predictive import perform_prediction_three_point_bending
import pytest
import numpy as np
import os
from pathlib import Path


def test_prediction():
    # check if KG is present
    kg_path = os.path.join(os.path.dirname(Path(__file__)),'prediction_data/prediction_kg.owl')
    if not os.path.exists(kg_path):
        raise Exception('The knowledge graph file is not present')

    pos_pred = perform_prediction_three_point_bending(knowledge_graph_file=kg_path, mode="cheap")

    # TODO carry units till the end!!!
    assert np.mean(pos_pred) == pytest.approx(125.00497593013911)
