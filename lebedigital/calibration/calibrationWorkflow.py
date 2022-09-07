# third party imports
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import rc

mpl.rcParams['font.size'] = 16
rc('text', usetex=True)

# local imports (probeye)
from probeye.definition.inverse_problem import InverseProblem
from probeye.definition.likelihood_model import GaussianLikelihoodModel
from probeye.ontology.knowledge_graph_export import export_knowledge_graph
from probeye.ontology.knowledge_graph_import import import_parameter_samples
from probeye.ontology.knowledge_graph_export import export_results_to_knowledge_graph
from probeye.inference.emcee.solver import EmceeSolver

# local imports (others)
from lebedigital.Prediction.three_point_bending_example import three_point_bending_example
from lebedigital.calibration.utils import query_KG, PosteriorPredictive
from lebedigital.calibration.forwardmodel_linear_elastic_cylinder import LinearElasticityCylinder


def perform_calibration(path_to_KG : str, calibrated_data_path : str, experiment_name='Wolf 8.2 Probe 1'):
    # =========================================
    #       Loading the experiment
    # =========================================
    # for concrete E 60 - 85 Gpa/60 - 85 10^9N/m2 /60-85 KN/mm^2
    #experiment_name = 'Wolf 8.2 Probe 1'
    # -- if KG is not there, using custom script
    # path_csv = '../MinimumWorkingExample/emodul/processed_data/Wolf 8.2 Probe 1.csv'
    # skip_last = 145
    # skip_init = 330
    #experiment_name = 'Wolf 8.2 Probe 1'
    # exp_output = load_experimental_data(experiment_name, skip_init, skip_last, path=path_csv)
    # -- if KG is present, comment the above
    #path_to_KG = '../usecases/MinimumWorkingExample/emodul/knowledge_graphs/'
    path = path_to_KG + '/' + experiment_name + '.ttl'
    exp_output = query_KG(path=path)
    plt.plot(exp_output['stress'], exp_output['displacement'] / exp_output['height'])  # checking loading data

    # This generated Latex errors!!! Maybe we require some extra package
    #plt.show()
    # ========================================
    #       Set Numerical values
    # ========================================
    # "uninformed" Normal prior for the E modulus
    loc_E = 100  # kN/mm2
    scale_E = 100  # this can be inferred too

    # Uniform prior for experimental noise/model discrepancy
    low_sigma = 0
    high_sigma = 0.005

    # =========================================
    #       Define the Inference Problem
    # =========================================

    # initialize the inverse problem with a useful name
    problem = InverseProblem("compression test calibration")

    # add all parameters to the problem
    problem.add_parameter('E', 'model',
                          tex="$E$",
                          info="Slope of the graph",
                          prior=('normal', {'mean': loc_E,
                                            'std': scale_E}))

    problem.add_parameter('sigma', 'likelihood',
                          tex=r"$\sigma$",
                          info="Std. dev, of 0-mean noise model",
                          prior=('uniform', {'low': low_sigma,
                                             'high': high_sigma}))

    # Load the forward model and add it to the Inverse problem
    linear_elasticity = LinearElasticityCylinder("linear_elasticity_cylinder")
    problem.add_forward_model(linear_elasticity)

    # ============================================
    #     Add test data to the Inverse Problem
    # ============================================

    y_test = exp_output['stress'] * (np.pi * (exp_output['diameter'] / 2) ** 2)  # in kN
    # add the experimental data
    problem.add_experiment(
        experiment_name,
        fwd_model_name="linear_elasticity_cylinder",
        sensor_values={
            'nu': 0.2,
            'height': exp_output['height'],
            'radius': exp_output['diameter'] / 2,
            'displacement_list': exp_output['displacement'],
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
    #dir_path = os.path.dirname(__file__)
    dir_path = calibrated_data_path
    basename_owl = os.path.basename(__file__).split(".")[0] + ".owl"
    knowledge_graph_file = os.path.join(dir_path, basename_owl)
    export_knowledge_graph(problem, knowledge_graph_file, data_dir=dir_path)

    # ========================================================================
    #        Solve problem with inference engine and write results to graph
    # ========================================================================

    # run inference step using emcee
    emcee_solver = EmceeSolver(problem, show_progress=True)
    inference_data = emcee_solver.run_mcmc(
        n_walkers=5,
        n_steps=40,
        n_initial_steps=5,
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
    E_pos = sample_dict['E'] * 1000  # N/mm2 ~ E_mean ~ 95E03N/mm2
    # three_point = three_point_bending_example()
    pos_pred = PosteriorPredictive(three_point_bending_example, known_input_solver=nu, parameter=E_pos)
    mean, sd = pos_pred.get_stats(samples=10)  # mean : ~365 N/mm2, sd = 30
    # ---- visualize posterior predictive
    posterior_pred_samples = pos_pred._samples
    # np.save('./post_pred.npy', posterior_pred_samples)
    plt.figure()
    sns.kdeplot(posterior_pred_samples)
    plt.xlabel('stress in x-direction')
    # more problems with latex when running doit
    #plt.show()
    #plt.savefig(calibrated_data_path + 'posterior_predictive.pdf',dpi =100)
    # plt.savefig('Figures/posterior_predictive.pdf', dpi=100)



# testing things
from pathlib import Path
#parent directory of the minimum working example
ParentDir = '../../usecases/MinimumWorkingExample'
exp_name = 'Wolf 8.2 Probe 1'
knowledge_graphs_directory = Path(ParentDir, 'emodul', 'knowledge_graphs')
calibrated_data_directory = Path(ParentDir, 'emodul', 'calibrated_data')
perform_calibration(str(knowledge_graphs_directory),calibrated_data_directory,exp_name)