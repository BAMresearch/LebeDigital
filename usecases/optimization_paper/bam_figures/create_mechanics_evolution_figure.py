import copy

import fenics_concrete
import matplotlib.pyplot as plt
import numpy as np
import pint
import pytest
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg

from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# plt.style.use('ggplot')
SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 20

rc('font', size=MEDIUM_SIZE)          # controls default text sizes
rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure titles


def create_mechanics_evolution_figure(input_parameter: dict, fig_path: str = "test_mechanics_evolution_plot.pdf"):
    #      add figure to tex file and snakemake workflow

    material_problem = fenics_concrete.ConcreteThermoMechanical()
    e_fkt = material_problem.mechanics_problem.E_fkt
    fc_fkt = material_problem.mechanics_problem.general_hydration_fkt

    # general parameters
    parameter = {
        "alpha_t": input_parameter["evoExAlphaT"],
        "alpha_tx": input_parameter["evoExAlphaTx"],
        "alpha_0": 0.0,
        "a_E": input_parameter["evoExaE"],
        "a_X": input_parameter["evoExafc"],
        "E": input_parameter["evoExE"],
        "X": input_parameter["evoExfc"],
    }

    alpha_list = np.arange(0, 1, 0.005)

    variation_dict = {
        "alpha_t": {
            "params": [parameter["alpha_t"], 0.0, 0.6],
            "fkt": e_fkt,
            "ylabel": "Elastic modulus $E$, GPa",
            "ylim": 60,
        },
        "a_E": {
            "params": [parameter["a_E"], 0.2, 1.3],
            "fkt": e_fkt,
            "ylabel": "Elastic modulus $E$, Gpa",
            "ylim": 60,
        },
        "a_X": {
            "params": [parameter["a_X"], 0.2, 1.3],
            "fkt": fc_fkt,
            "ylabel": "Compressive strength $f_c$, MPa",
            "ylim": 40,
        },
    }

    # setup plot
    fig, axs = plt.subplots(1, len(variation_dict), figsize=(5 * len(variation_dict), 4))
    ureg.setup_matplotlib()

    i = 0
    temp_key = {"alpha_t": r"$\alpha_t$", "a_E": r"$a_E$", "a_X": r"$a_X$"}
    for key in variation_dict.keys():
        p = copy.deepcopy(parameter)
        var_par_list = variation_dict[key]["params"]
        fkt = variation_dict[key]["fkt"]

        for value in var_par_list:
            p[key] = value
            y_list = []
            for alpha in alpha_list:
                y_list.append(fkt(alpha, p))
            #axs[i].plot(alpha_list, y_list, label=key + " = " + str(value))
            axs[i].plot(alpha_list, y_list, label=temp_key[key] + " = " + str(value))
            axs[i].legend()
            axs[i].set_xlabel(f"Degree of hydration $\\alpha$")
            axs[i].set_ylabel(variation_dict[key]["ylabel"])
            axs[i].set_xlim([0, 1])
            axs[i].set_ylim([0, variation_dict[key]["ylim"]])

        i += 1

    fig.tight_layout()
    # plt.show()
    fig.savefig(fig_path)


if __name__ == "__main__":
    parameter = {
        "evoExAlphaT": 0.2,
        "evoExAlphaTx": 0.8,
        "evoExaE": 0.5,
        "evoExafc": 0.5,
        "evoExE": 50,
        "evoExfc": 30,
    }

    create_mechanics_evolution_figure(parameter)