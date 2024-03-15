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


def create_heat_release_figure(parameter: dict, fig_path: str = "test_heat_realease_plot.pdf"):
    #      add figure to tex file and snakemake workflow

    parameter["B1"] = parameter["heatExBOne"]
    parameter["B2"] = parameter["heatExBTwo"]
    parameter["eta"] = parameter["heatExEta"]
    parameter["alpha_max"] = parameter["heatExAlphaMax"]
    parameter["E_act"] = parameter["heatExEAct"]
    parameter["T_ref"] = parameter["heatExTRef"]
    parameter["Q_pot"] = parameter["heatExQPot"]

    T = parameter["heatExT"]  # temperature...
    dt = parameter["heatExDt"]  # dt
    time_total = parameter["heatExTimeTotal"]
    # what does T and dt do???

    time_list = np.arange(0, time_total, dt)

    variation_dict = {
        "B1": [parameter["B1"], 2.0e-4, 3.7e-4],
        "B2": [parameter["B2"], 0.0001, 0.01],
        "eta": [parameter["eta"], 9, 4.5],
        "Q_pot": [parameter["Q_pot"], 350e3, 650e3],
    }

    material_problem = fenics_concrete.ConcreteThermoMechanical()
    hydration_fkt = material_problem.get_heat_of_hydration_ftk()

    fig, axs = plt.subplots(2, len(variation_dict), figsize=(20, 7))
    ureg.setup_matplotlib()

    i = 0
    for key in variation_dict.keys():
        p = copy.deepcopy(parameter)
        for value in variation_dict[key]:
            p[key] = value
            heat_list, doh_list = hydration_fkt(T, time_list, dt, p)

            delta_heat = np.diff(heat_list) / dt
            plot_time = time_list[:-1]
            # add pint units to plot_time and delta_heat
            plot_time = plot_time * ureg.second
            # time_list = time_list * ureg.second
            delta_heat = delta_heat * ureg.watt / ureg.kg
            heat_list = heat_list * ureg.joule / ureg.kg

            # convert plot_time to hours
            plot_time = plot_time.to(ureg.hour)
            # time_list = time_list.to(ureg.hour)
            # convert delta_heat to mW/kg
            delta_heat = delta_heat.to(ureg.mW / ureg.kg)
            # heat_list = heat_list.to(ureg.mW / ureg.kg)

            axs[0][i].set_ylim([0, 6])
            axs[0][i].set_xlim([0, 24])
            # plot delta heat over time with a legend
            axs[0][i].plot(plot_time, delta_heat, label=key + " = " + str(value))

            # cummulative heat release
            axs[1][i].set_ylim([0, 400])
            axs[1][i].set_xlim([0, 24 * 4])
            axs[1][i].plot(plot_time, heat_list[:-1], label=key + " = " + str(value))

            # plt.plot(time_list, heat_list, label=parameter + " = " + str(value))
        axs[0][i].legend()
        # set legend to lower right corner
        axs[1][i].legend()
        axs[1][i].legend(loc="lower right")

        #axs[0][i].set_ylabel(f"Heat release rate in {delta_heat.units}")
        axs[0][i].set_ylabel(f"Heat release rate $(mW/kg)$")
        axs[0][i].set_xlabel(f"time in {plot_time.units}")
        #axs[1][i].set_ylabel(f"Cumulated heat release in {heat_list.units}")
        axs[1][i].set_ylabel(f"Cumulated heat release $(j/kg)$")
        axs[1][i].set_xlabel(f"time in {plot_time.units}")

        i += 1

    fig.tight_layout()
    # plt.show()

    fig.savefig(fig_path)


if __name__ == "__main__":
    parameter = {
        "heatExBOne": 3e-4,
        "heatExBTwo": 0.001,
        "heatExEta": 6,
        "heatExAlphaMax": 0.875,
        "heatExEAct": 47002,
        "heatExTRef": 25,
        "heatExQPot": 500e3,
        "heatExT": 20,
        "heatExDt": 60,
        "heatExTimeTotal": 60 * 60 * 24 * 4,
    }

    create_heat_release_figure(parameter)