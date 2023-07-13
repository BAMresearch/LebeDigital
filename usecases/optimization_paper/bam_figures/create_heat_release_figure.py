import copy

import fenics_concrete
import matplotlib.pyplot as plt
import numpy as np
import pint
import pytest
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg


def create_heat_release_figure():
    # todo add units for the input, convert to correct untis
    #      add units for the output
    #      change time to hours our days in plot
    #      add legend
    #      add title
    #      add axis labels
    #      add figure to tex file and snakemake workflow

    T = 20  # temperature...
    dt = 60  # dt
    time_total = 100000
    # what does T and dt do???

    time_list = np.arange(0, time_total, dt)
    parameter = {}
    parameter["B1"] = 3e-4
    parameter["B2"] = 0.001
    parameter["eta"] = 6
    parameter["alpha_max"] = 0.875
    parameter["E_act"] = 47002
    parameter["T_ref"] = 25
    parameter["Q_pot"] = 500e3

    variation_dict = {
        "B1": [parameter["B1"], 2.0e-4, 3.7e-4],
        "B2": [parameter["B2"], 0.0001, 0.01],
        "eta": [parameter["eta"], 9, 4.5],
    }

    material_problem = fenics_concrete.ConcreteThermoMechanical()
    hydration_fkt = material_problem.get_heat_of_hydration_ftk()

    fig, axs = plt.subplots(1, len(variation_dict), figsize=(15, 3))
    ureg.setup_matplotlib()

    i = 0
    for key in variation_dict.keys():
        p = copy.deepcopy(parameter)
        for value in variation_dict[key]:
            p[key] = value
            heat_list, doh_list = hydration_fkt(T, time_list, dt, p)

            delta_heat = np.diff(heat_list) / dt
            plot_time = time_list[:-1]

            axs[i].set_ylim([0, 0.006])
            axs[i].plot(plot_time, delta_heat)

            # plt.plot(time_list, heat_list, label=parameter + " = " + str(value))

        i += 1

    # initiate material problem
    # get the respective function
    #
    # heat_list, doh_list = hydration_fkt(T, time_list, dt, parameter)
    # print(heat_list)
    # print(doh_list)
    # print(time_list)

    # #ureg.setup_matplotlib(enable=False)
    # ureg2 = pint.UnitRegistry(auto_reduce_dimensions=True)

    # fig.suptitle('Influence of aggregate ratio on effective concrete properties')
    # axs[0].plot(time_list, heat_list)
    # # axs[0].set_ylabel(f"Young's modulus in {E_list.units}")
    # # axs[0].set_xlabel('aggregate volume fraction')
    # axs[1].plot(time_list, doh_list)
    # # axs[1].set_ylabel(f"Poission's ratio")
    # # axs[1].set_xlabel('aggregate volume fraction')
    # axs[2].plot(time_list, delta_heat)
    # axs[2].set_ylabel(f"compressive strength {fc_list.units}")
    # axs[2].set_xlabel('aggregate volume fraction')

    # plt.subplots_adjust(wspace=0.4)
    fig.tight_layout()
    plt.show()

    # fig.savefig(fig_path)


if __name__ == "__main__":
    create_heat_release_figure()
