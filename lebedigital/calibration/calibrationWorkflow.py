# third party imports
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import rc

mpl.rcParams["font.size"] = 16
rc("text", usetex=True)

# local imports (probeye)
from probeye.definition.inverse_problem import InverseProblem
from probeye.definition.likelihood_model import GaussianLikelihoodModel
from probeye.ontology.knowledge_graph_export import export_knowledge_graph
from probeye.ontology.knowledge_graph_import import import_parameter_samples
from probeye.ontology.knowledge_graph_export import export_results_to_knowledge_graph
from probeye.inference.emcee.solver import EmceeSolver

# local imports (others)
from lebedigital.calibration.forwardmodel_linear_elastic_cylinder import (
    LinearElasticityCylinder,
)


def _test_E_mod_calibration_metadata(calibration_metadata: dict):
    """
    Checks the calibration metadata dictionary for completeness
    Parameters
    ----------
    calibration_metadata :

    Returns
    -------

    """
    if {"E_loc", "E_scale"} <= calibration_metadata.keys():
        return True
    else:
        return False


def _test_E_mod_experimental_data(experimental_data: dict):
    """
    Performs a sanity check for the experimental data dict.
    Parameters
    ----------
    experimental_data :

    Returns
    -------

    """

    if {"exp_name", "force", "displacement", "height", "diameter"} <= experimental_data.keys():
        return True
    else:
        return False


def esimate_Youngs_modulus(
        experimental_data: dict, calibration_metadata: dict, calibrated_data_path: str, mode="full"
):
    """
    Function to solve an inverse problem using Bayesian inference to infer Young's Modulus (E), with experimental
    load-displacement being the observed data and FE based concrete compression solver.
    The inferred E modulus is used in a posterior predictive setting, to propagate the uncertainity and compute the
    stress in 3 point bending test.

    Parameters
    ----------
    experimental_data : dict
        Must contain the following keys:
        - exp_name (str) : The name of the experiment
        - force : The force value (kN) array
        - displacement : The displacement array (mm)
        - height: The height of the sample (mm)
        - diameter : The diameter of the sample (mm)
    calibration_metadata : dict
        Parameters needed to perform the calibration is passed here. Note that this will work
        for this specific setup for E modulus calibration. Following keys must be present.
        - E_loc : The mean value of the E mod for the prior (initial guess)
        - E_scale : The s.d value of the E mod for the prior (initial guess)
    calibrated_data_path : str
        Path where the calibrated results needs to be stored. The calibration results along with the inverse problem
        setting is stored in this path as knowledge graph
    mode : "full" or "cheap". For testing purposes.

    Returns
    -------
    E_pos : np.array
        The samples of the inferred Young's Modulus.
    posterior_pred_samples: np.array
        The posterior predictive values of the stress in the three point bending beam obtained using the inferred E-modulus.

    """

    # =========================================
    #       Loading the experiment
    # =========================================
    # perform check of the experimental data keys
    assert (
            _test_E_mod_experimental_data(experimental_data) == True
    ), "Some values are missing in the experimental data"
    exp_output = experimental_data  #

    # ========================================
    #       Set Numerical values
    # ========================================
    # performing check of calibration metadata
    assert (
            _test_E_mod_calibration_metadata(calibration_metadata) == True
    ), "Some values are missing in the calibration metadata"
    # "uninformed" Normal prior for the E modulus
    loc_E = calibration_metadata["E_loc"]  # 100  kN/mm2
    scale_E = calibration_metadata["E_scale"]  # 100  this can be inferred too

    # Uniform prior for experimental noise/model discrepancy
    low_sigma = 0
    high_sigma = 0.005

    # =========================================
    #       Define the Inference Problem
    # =========================================

    # initialize the inverse problem with a useful name
    problem = InverseProblem("compression test calibration")

    # add all parameters to the problem
    problem.add_parameter(
        "E",
        "model",
        tex="$E$",
        info="Slope of the graph",
        prior=("normal", {"mean": loc_E, "std": scale_E}),  # can add log_normal here
    )

    problem.add_parameter(
        "sigma",
        "likelihood",
        tex=r"$\sigma$",
        info="Std. dev, of 0-mean noise model",
        prior=("uniform", {"low": low_sigma, "high": high_sigma}),
    )

    # Load the forward model and add it to the Inverse problem
    linear_elasticity = LinearElasticityCylinder("linear_elasticity_cylinder")
    problem.add_forward_model(linear_elasticity)

    # ============================================
    #     Add experimental data to the Inverse Problem
    # ============================================
    experiment_name = exp_output["exp_name"]
    y_test = exp_output["force"]  # in kN
    # add the experimental data
    problem.add_experiment(
        experiment_name,
        fwd_model_name="linear_elasticity_cylinder",
        sensor_values={
            "nu": 0.2,
            "height": exp_output["height"],
            "radius": exp_output["diameter"] / 2,
            "displacement_list": exp_output["displacement"],
            linear_elasticity.output_sensor.name: y_test,
        },
    )

    # ==============================================
    #       Add likelihood model(s)
    # ==============================================

    # add the likelihood model to the problem
    problem.add_likelihood_model(
        GaussianLikelihoodModel(
            prms_def="sigma",
            experiment_name=experiment_name,
            model_error="additive",
            name="SimpleLikelihoodModel",
        )
    )

    # give problem overview
    problem.info()

    # ==============================================
    #        Export knowledge graph
    # ==============================================

    # create the knowledge graph and print it to file
    # dir_path = os.path.dirname(__file__)
    dir_path = calibrated_data_path
    # basename_owl = os.path.basename(__file__).split(".")[0] + ".owl"
    basename_owl = os.path.basename(__file__).split(".")[0] + exp_output['exp_name'] + ".owl"
    knowledge_graph_file = os.path.join(dir_path, basename_owl)
    export_knowledge_graph(problem, knowledge_graph_file, data_dir=dir_path)

    # ========================================================================
    #        Solve problem with inference engine and write results to graph
    # ========================================================================

    # run inference step using emcee
    emcee_solver = EmceeSolver(problem, seed=0, show_progress=True)

    if mode == "cheap":
        inference_data = emcee_solver.run_mcmc(
            n_walkers=3,
            n_steps=4,
            n_initial_steps=2,
        )
    else:
        inference_data = emcee_solver.run_mcmc(
            n_walkers=5,
            n_steps=50,
            n_initial_steps=10,
        )

    # export the results from the 'inference_data' object to the graph
    export_results_to_knowledge_graph(
        problem,
        inference_data,
        knowledge_graph_file,
        data_dir=dir_path,
    )

    # =======================================================================
    #                          Query the knowledge graph
    # =======================================================================

    # get the samples from the knowledge graph
    sample_dict = import_parameter_samples(knowledge_graph_file)
    E_pos = sample_dict["E"]  # kN/mm2 ~ E_mean ~ 95E03N/mm2

    return E_pos
