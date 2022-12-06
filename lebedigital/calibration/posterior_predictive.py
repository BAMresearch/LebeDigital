# third party imports
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import rc

# local imports (probeye)
from probeye.ontology.knowledge_graph_import import import_parameter_samples


# local imports (others)
from lebedigital.Prediction.three_point_bending_example import (
    three_point_bending_example,
)
from lebedigital.calibration.utils import (
    PosteriorPredictive,
)

def perform_prediction_three_point_bending(knowledge_graph_file:str,forward_solver=three_point_bending_example,mode = 'cheap'):
    """

    Parameters
    ----------
    knowledge_graph_file : (str) A KG file containing the stats of the inferred parameter.
    forward_solver : The solver through which parametric uncertainty needs to be propagated
    mode : "full" or "cheap". For testing purposes.

    Returns
    -------
    posterior_pred_samples: np.array (N/mm2)
    The posterior predictive values of the stress in the three point bending beam obtained using the inferred E-modulus.
    """


    # =======================================================================
    #                          Query the knowledge graph
    # =======================================================================

    # get the samples from the knowledge graph
    sample_dict = import_parameter_samples(knowledge_graph_file)

    # ========================================================================
    #       Posterior Predictive
    # ========================================================================
    nu = 0.2
    # conversion from kN/mm^2 to N/mm2
    E_pos = sample_dict["E"] * 1000  # N/mm2 ~ E_mean ~ 95E03N/mm2, currently 'E' and the unit conversion factor 1000
    # is harcoded here.
    # three_point = three_point_bending_example()
    pos_pred = PosteriorPredictive(
        forward_solver, known_input_solver=nu, parameter=E_pos
    )

    if mode == "cheap":
        mean, sd = pos_pred.get_stats(samples=1)  # mean : ~365 N/mm2, sd = 30
    else:
        mean, sd = pos_pred.get_stats(samples=10)  # mean : ~365 N/mm2, sd = 30
    # ---- visualize posterior predictive
    posterior_pred_samples = pos_pred._samples
    # np.save('./post_pred.npy', posterior_pred_samples)
    plt.figure()
    sns.kdeplot(posterior_pred_samples)
    plt.xlabel("stress in x-direction")

    return posterior_pred_samples

