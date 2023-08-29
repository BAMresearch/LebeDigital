import matplotlib.pyplot as plt
import numpy as np

from lebedigital.demonstrator_scripts.beam_design import check_beam_design
from lebedigital.unit_registry import ureg


def simple_setup(input_parameter, height, fc, load):
    """
    testing some


    Parameters
    ----------
    input_parameter : dict
        dictionary containing the constant input parameters
    height : pint quantity
        height of the beam
    fc : pint quantity
        compressive strength of the concrete
    load : pint quantity
        load on the beam

    Returns
    -------
    out : dict
        dictionary containing the output parameters
    """

    out = check_beam_design(
        span=input_parameter["beamExSpan"] * ureg(input_parameter["beamExSpanUnit"]),
        width=input_parameter["beamExWidth"] * ureg(input_parameter["beamExWidthUnit"]),
        height=height,
        point_load=load,
        distributed_load=0 * ureg("N/mm"),
        compr_str_concrete=fc,
        yield_str_steel=input_parameter["beamExYieldStrSteel"] * ureg(input_parameter["beamExYieldStrSteelUnit"]),
        steel_dia_bu=input_parameter["beamExSteelDiaBu"] * ureg(input_parameter["beamExSteelDiaBuUnit"]),
        cover_min=input_parameter["beamExCoverMin"] * ureg(input_parameter["beamExCoverMinUnit"]),
    )

    return out


def beam_design_plot(input_parameter, n, fig_path: str = "beam_design_plot.pdf", display_output: bool = False):
    min_fc = 1
    max_fc = 100
    fc_unit = "N/mm^2"
    min_load = 10
    max_load = 200
    load_unit = "kN"
    min_height = 150
    max_height = 500
    height_unit = "mm"
    crosssection_unit = "cm^2"

    def get_list(min, max, n):
        step = (max - min) / (n - 1)
        return np.arange(min, max + step, step)

    def get_plot_lists(order_list, values_dict, input_parameter):
        if display_output:
            print(f"compute: {order_list[0]} - {order_list[1]}")

        x_list = values_dict[order_list[0]]["list"]
        y_list = values_dict[order_list[1]]["list"]
        constant = values_dict[order_list[2]]["constant"]

        order_list = [0 if item == "height" else item for item in order_list]
        order_list = [1 if item == "fc" else item for item in order_list]
        order_list = [2 if item == "load" else item for item in order_list]

        crosssections = np.zeros(shape=(len(x_list), len(y_list)))
        fc_errors = np.zeros(shape=(len(x_list), len(y_list)))
        A_errors = np.zeros(shape=(len(x_list), len(y_list)))
        constraint = np.zeros(shape=(len(x_list), len(y_list)))

        values = [0, 0, 0]
        values[order_list[2]] = constant
        max_admissable_area = 0.0
        for i, x in enumerate(x_list):
            values[order_list[0]] = x
            for j, y in enumerate(y_list):
                values[order_list[1]] = y

                out = simple_setup(input_parameter, *values)
                crosssection = out["crosssection"]
                crosssection.ito(crosssection_unit)
                crosssection = crosssection.magnitude
                crosssections[i][j] = crosssection
                fc_errors[i][j] = out["constraint_min_fc"]
                A_errors[i][j] = out["constraint_max_steel_area"]
                constraint[i][j] = out["constraint_beam_design"]

                if constraint[i][j] <= 0:
                    if max_admissable_area < crosssections[i][j]:
                        max_admissable_area = crosssections[i][j]

        return crosssections, fc_errors, A_errors, constraint, max_admissable_area

    def plot_contour(ax, x, y, Z, color, linewidth, linestyle="solid"):
        max = np.nanmax(Z)
        min = np.nanmin(Z)

        levels = [min - 1, 0.0, max + 1]

        X, Y = np.meshgrid(y, x)
        ax.contour(
            X.magnitude,
            Y.magnitude,
            Z,
            colors=color,
            levels=levels,
            linewidths=linewidth,
            linestyles=linestyle,
        )

    def plot_contour_filled_white(ax, x, y, Z):
        max = np.nanmax(Z)
        min = np.nanmin(Z)
        levels = [0.0, max]

        X, Y = np.meshgrid(y, x)
        ax.contourf(X.magnitude, Y.magnitude, Z, colors="white", levels=levels)

    def plot_contour_filled(ax, x, y, Z, xlabel, ylabel, title, colorbar=False, max=None, min=None):
        if not max:
            max = np.nanmax(Z)
        if not min:
            min = 0
        max = max * 1.1

        n_levels = 8
        level_step = max / (n_levels - 1)
        levels = np.arange(min, max + level_step, level_step)

        X, Y = np.meshgrid(y, x)
        cs = ax.contourf(X.magnitude, Y.magnitude, Z, vmin=min, vmax=max - level_step, levels=levels)
        ax.set_ylabel(xlabel)
        ax.set_xlabel(ylabel)
        ax.set_title(title)

        if colorbar:
            cbar = plt.colorbar(cs)
            cbar.set_label(f"total steel cross section in {ureg(crosssection_unit).units}")

    values = {
        "fc": {
            "list": get_list(min_fc, max_fc, n) * ureg(fc_unit),
            "constant": input_parameter["beamExComprStrConcreteC"]
            * ureg(input_parameter["beamExComprStrConcreteCUnit"]),
        },
        "load": {
            "list": get_list(min_load, max_load, n) * ureg(load_unit),
            "constant": input_parameter["beamExPointLoadC"] * ureg(input_parameter["beamExPointLoadCUnit"]),
        },
        "height": {
            "list": get_list(min_height, max_height, n) * ureg(height_unit),
            "constant": input_parameter["beamExHeightC"] * ureg(input_parameter["beamExHeightCUnit"]),
        },
    }

    chosen_plots = [
        ("fc", "load", "height"),
        ("height", "load", "fc"),
        ("height", "fc", "load"),
    ]

    # empty lists
    crosssection_plots = [0] * len(chosen_plots)
    fc_error_plots = [0] * len(chosen_plots)
    A_error_plots = [0] * len(chosen_plots)
    constraint_plots = [0] * len(chosen_plots)
    plot_values = {}
    max_plot_area = 0.0

    # loop over all three plots
    for i, plot in enumerate(chosen_plots):
        (
            crosssection_plots[i],
            fc_error_plots[i],
            A_error_plots[i],
            constraint_plots[i],
            max_admissible_area,
        ) = get_plot_lists(plot, values, input_parameter)
        if max_admissible_area > max_plot_area:
            max_plot_area = max_admissible_area

        plot_values[plot] = {
            "crosssections": crosssection_plots[i],
            "fc_errors": fc_error_plots[i],
            "A_errors": A_error_plots[i],
            "constraint": constraint_plots[i],
        }

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    max_crosssection = 0.0
    for plot in chosen_plots:
        plot_max = np.nanmax(plot_values[plot]["crosssections"])
        if plot_max > max_crosssection:
            max_crosssection = plot_max

    max_crosssection = max_plot_area

    for i, plot in enumerate(chosen_plots):
        colorbar_bool = False
        if i + 1 == len(chosen_plots):
            colorbar_bool = True

        plot_contour_filled(
            axs[i],
            values[plot[0]]["list"],
            values[plot[1]]["list"],
            plot_values[plot]["crosssections"],
            xlabel=f'{plot[0]} in {values[plot[0]]["list"].units}',
            ylabel=f'{plot[1]} in {values[plot[1]]["list"].units}',
            title=f'{plot[2]} = {values[plot[2]]["constant"]}',
            colorbar=colorbar_bool,
            max=max_crosssection,
        )

        plot_contour_filled_white(
            axs[i],
            values[plot[0]]["list"],
            values[plot[1]]["list"],
            plot_values[plot]["constraint"],
        )

        plot_contour(
            axs[i],
            values[plot[0]]["list"],
            values[plot[1]]["list"],
            plot_values[plot]["fc_errors"],
            color="black",
            linewidth=3.0,
            linestyle="dashed",
        )

        plot_contour(
            axs[i],
            values[plot[0]]["list"],
            values[plot[1]]["list"],
            plot_values[plot]["A_errors"],
            color="black",
            linewidth=3.0,
            linestyle="dotted",
        )

    fig.tight_layout()
    fig.savefig(fig_path)


if __name__ == "__main__":
    beam_design_example_parameters = {
        "beamExSpan": 675,
        "beamExSpanUnit": "cm",
        "beamExWidth": 200,
        "beamExWidthUnit": "mm",
        "beamExYieldStrSteel": 500,
        "beamExYieldStrSteelUnit": "N/mm^2",
        "beamExSteelDiaBu": 10,
        "beamExSteelDiaBuUnit": "mm",
        "beamExCoverMin": 2.5,
        "beamExCoverMinUnit": "cm",
        "beamExHeightC": 450,
        "beamExHeightCUnit": "mm",
        "beamExPointLoadC": 50,
        "beamExPointLoadCUnit": "kN",
        "beamExComprStrConcreteC": 40,
        "beamExComprStrConcreteCUnit": "N/mm^2",
    }

    beam_design_plot(beam_design_example_parameters, n=10, display_output=True)
