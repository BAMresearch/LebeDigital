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
from lebedigital.Prediction.three_point_bending_example import (
    three_point_bending_example,
)
from lebedigital.calibration.utils import (
    PosteriorPredictive,
)
from lebedigital.calibration.forwardmodel_linear_elastic_cylinder import (
    LinearElasticityCylinder,
)


def _test_E_mod_calibration_metadata(calibration_metadata: dict):
    """
    Checks the calibration metadata dictionary for completeness
    Args:
        calibration_metadata ():

    Returns:

    """
    if {"E_loc", "E_scale"} <= calibration_metadata.keys():
        return True
    else:
        return False


def _test_E_mod_experimental_data(experimental_data: dict):
    """
    Performs a sanity check for the experimental data dict.
    Args:
        experimental_data ():

    Returns:

    """
    if {"exp_name", "force", "displacement", "height", "diameter"} <= experimental_data.keys():
        return True
    else:
        return False


def esimate_Youngs_modulus(
        experimental_data: dict, calibration_metadata: dict, calibrated_data_path: str, mode="full"
):
    """

    Args:
        calibrated_data_path (): Path where the calibrated results needs to be stored
        experimental_data (): Must contain the keys "exp_name" (str),"force"(kN), "displacement","height"(mm) and "diameter"(mm). Force
        and displacement needs to be arrays (The data stored in KG should ensure that it is from the third loading cycle).
         These needs to be computed externally and provided here by the knowledge graph module.
        calibration_metadata (): Parameters needed to perform the calibration is passed here. Note that this will work
        for this specific setup for E modulus calibration.
        mode (): "full" or "cheap". cheap is used just for tests

    Returns:
        E_pos : Posterior samples of the E mod.
        posterior_pred_samples : Posterior predictive samples
        The script also saves a Knowledge graph at the path specified by the "calibration_data_path"

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
    #     Add test data to the Inverse Problem
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
    basename_owl = os.path.basename(__file__).split(".")[0] + exp_output['exp_name']
    knowledge_graph_file = os.path.join(dir_path, basename_owl)
    export_knowledge_graph(problem, knowledge_graph_file, data_dir=dir_path)

    # ========================================================================
    #        Solve problem with inference engine and write results to graph
    # ========================================================================

    # run inference step using emcee
    emcee_solver = EmceeSolver(problem, seed=0, show_progress=True)

    if mode == "cheap":
        inference_data = emcee_solver.run_mcmc(
            n_walkers=4,
            n_steps=1,
            n_initial_steps=4,
        )
    else:
        inference_data = emcee_solver.run_mcmc(
            n_walkers=5,
            n_steps=100,
            n_initial_steps=20,
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

    # ========================================================================
    #       Posterior Predictive
    # ========================================================================
    nu = 0.2
    E_pos = sample_dict["E"] * 1000  # N/mm2 ~ E_mean ~ 95E03N/mm2
    # three_point = three_point_bending_example()
    pos_pred = PosteriorPredictive(
        three_point_bending_example, known_input_solver=nu, parameter=E_pos
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

    return E_pos, posterior_pred_samples
