import matplotlib.pyplot as plt
import numpy as np
import pint
import pytest
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.simulation.concrete_homogenization import concrete_homogenization
from lebedigital.unit_registry import ureg


def create_homogenization_figure(parameter_with_units: dict, fig_path: str = "homogenization_plot.pdf"):
    # initialize dictionary
    parameters = {}

    # paste data
    parameters["paste_E"] = float(parameter_with_units["pasteE"]) * ureg(parameter_with_units["pasteEunit"])
    parameters["paste_nu"] = float(parameter_with_units["pastenu"]) * ureg(parameter_with_units["pastenuunit"])
    parameters["paste_C"] = float(parameter_with_units["pasteC"]) * ureg(parameter_with_units["pasteCunit"])
    parameters["paste_kappa"] = float(parameter_with_units["pastekappa"]) * ureg(
        parameter_with_units["pastekappaunit"]
    )
    parameters["paste_rho"] = float(parameter_with_units["pasterho"]) * ureg(parameter_with_units["pasterhounit"])
    parameters["paste_fc"] = float(parameter_with_units["pastefc"]) * ureg(parameter_with_units["pastefcunit"])
    parameters["paste_Q"] = float(parameter_with_units["pasteQ"]) * ureg(parameter_with_units["pasteQunit"])

    # aggregate data
    parameters["aggregates_E"] = float(parameter_with_units["aggregatesE"]) * ureg(
        parameter_with_units["aggregatesEunit"]
    )
    parameters["aggregates_nu"] = float(parameter_with_units["aggregatesnu"]) * ureg(
        parameter_with_units["aggregatesnuunit"]
    )
    parameters["aggregates_C"] = float(parameter_with_units["aggregatesC"]) * ureg(
        parameter_with_units["aggregatesCunit"]
    )
    parameters["aggregates_kappa"] = float(parameter_with_units["aggregateskappa"]) * ureg(
        parameter_with_units["aggregateskappaunit"]
    )
    parameters["aggregates_rho"] = float(parameter_with_units["aggregatesrho"]) * ureg(
        parameter_with_units["aggregatesrhounit"]
    )

    # breakpoint()

    for key in parameters.keys():
        assert parameters[key] == parameters[key]

    n_points = 100
    step_size = 1 / n_points
    aggregates_vol_frac_list = np.arange(0, 1 + step_size, step_size)
    assert aggregates_vol_frac_list.max() <= 1.0

    E_list = []
    nu_list = []
    fc_list = []
    kappa_list = []
    Q_list = []

    for value in aggregates_vol_frac_list:
        parameters["aggregates_vol_frac"] = value * ureg("dimensionless")

        results = concrete_homogenization(parameters)

        E_list.append(results["E"].magnitude)
        nu_list.append(results["nu"].magnitude)
        fc_list.append(results["fc"].magnitude)
        kappa_list.append(results["kappa"].magnitude)
        Q_list.append(results["Q"].magnitude)

    E_list = E_list * results["E"].units
    nu_list = nu_list * results["nu"].units
    fc_list = fc_list * results["fc"].units
    kappa_list = kappa_list * results["kappa"].units
    Q_list = Q_list * results["Q"].units

    # #ureg.setup_matplotlib(enable=False)
    # ureg2 = pint.UnitRegistry(auto_reduce_dimensions=True)
    ureg.setup_matplotlib()

    fig, axs = plt.subplots(1, 5, figsize=(15, 3))
    # fig.suptitle('Influence of aggregate ratio on effective concrete properties')
    axs[0].plot(aggregates_vol_frac_list, E_list)
    axs[0].set_ylabel(f"Young's modulus in {E_list.units}")
    axs[0].set_xlabel("aggregate volume fraction")
    axs[1].plot(aggregates_vol_frac_list, nu_list)
    axs[1].set_ylabel(f"Poission's ratio")
    axs[1].set_xlabel("aggregate volume fraction")
    axs[2].plot(aggregates_vol_frac_list, fc_list)
    axs[2].set_ylabel(f"compressive strength {fc_list.units}")
    axs[2].set_xlabel("aggregate volume fraction")
    axs[3].plot(aggregates_vol_frac_list, kappa_list)
    axs[3].set_ylabel(f"thermal conductivity\nin {kappa_list.units}")
    axs[3].set_xlabel("aggregate volume fraction")
    axs[4].plot(aggregates_vol_frac_list, Q_list)
    axs[4].set_ylabel(f"Heat release {Q_list.units}")
    axs[4].set_xlabel("aggregate volume fraction")

    # plt.subplots_adjust(wspace=0.4)
    fig.tight_layout()

    fig.savefig(fig_path)


if __name__ == "__main__":
    homogenization_example_parameters = {
        "pasteE": 30e9,
        "pasteEunit": "Pa",
        "pastenu": 0.2,
        "pastenuunit": "dimensionless",
        "pasteC": 870,
        "pasteCunit": "J/kg/K",
        "pastekappa": 1.8,
        "pastekappaunit": "W/m/K",
        "pasterho": 2400,
        "pasterhounit": "kg/m^3",
        "pastefc": 30e6,
        "pastefc_unit": "Pa",
        "pasteQ": 250000,
        "pasteQunit": "J/kg",
        "aggregatesE": 25e9,
        "aggregatesE_unit": "Pa",
        "aggregatesnu": 0.3,
        "aggregatesnu_unit": "dimensionless",
        "aggregatesC": 840,
        "aggregatesCunit": "J/kg/K",
        "aggregateskappa": 0.8,
        "aggregateskappaunit": "W/m/K",
        "aggregatesrho": 2600,
        "aggregatesrhounit": "kg/m^3",
    }

    create_homogenization_figure(homogenization_example_parameters)
