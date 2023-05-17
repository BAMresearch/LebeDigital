# third party imports
import os
from pathlib import Path
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import rc
import fenics_concrete

# local imports (probeye)
from probeye.ontology.knowledge_graph_import import import_parameter_samples

# local imports (others)
from lebedigital.simulation.three_pont_bending_beam import three_point_bending_example
from lebedigital.calibration.utils import PosteriorPredictive
from lebedigital.unit_registry import ureg


def wrapper_three_point_bending(parameter, known_input):
    """
    This is standard template any forward solver needs to be wrapped in
    Parameters
    ----------
    parameter : parameter is E in kN/mm2
    known_input : Poisson's ratio 'nu'

    Returns
    -------
    y : magnitude of stress output from the three point bending test

    """
    parameters = fenics_concrete.Parameters()
    parameters['E'] = parameter * 1000 * ureg('N/mm^2')
    parameters['nu'] = known_input * ureg('')

    y = three_point_bending_example(parameters)
    return y.magnitude


def perform_prediction(forward_solver:callable, parameter:list, nu:float = 0.2, no_sample :int = 50, mode ='cheap', viz=False):
    """

    Parameters
    ----------
    forward_solver : (callable) The solver through which parametric uncertainty needs to be propagated
    parameter : (list) the samples of the parameter which was calibrated. (E for eg)
    nu: Known input to the solver ie. nu
    no_sample: total no of samples for the MC estimate
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
    # sample_dict = import_parameter_samples(knowledge_graph_file)'
    # ignoring the knowledge graph

    # ========================================================================
    #       Posterior Predictive
    # ========================================================================
    # nu = known_input
    # E_pos = np.array(parameter)   # N/mm2 ~ E_mean ~ 95E03N/mm2, currently 'E' and the unit conversion factor 1000
    # is harcoded here.
    pos_pred = PosteriorPredictive(
        forward_solver, known_input_solver=nu, parameter=np.array(parameter)
    )

    if mode == "cheap":
        mean, sd = pos_pred.get_stats(samples=5)  # mean : ~365 N/mm2, sd = 30
    else:
        mean, sd = pos_pred.get_stats(samples=no_sample)  # mean : ~365 N/mm2, sd = 30
    # ---- visualize posterior predictive
    posterior_pred_samples = pos_pred._samples
    # np.save('./post_pred.npy', posterior_pred_samples)
    if viz:
        plt.figure()
        sns.kdeplot(posterior_pred_samples)
        plt.xlabel("stress in x-direction")

    return posterior_pred_samples

if __name__ == "__main__":
    E_samples = [30.1,30.2,30.53,30.8,29.6]
    random.shuffle(E_samples)
    pos_pred = perform_prediction(forward_solver=wrapper_three_point_bending,
                                  parameter=E_samples)
    print(pos_pred)