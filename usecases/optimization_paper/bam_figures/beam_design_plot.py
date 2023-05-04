from lebedigital.demonstrator_scripts.beam_design import check_beam_design
from lebedigital.unit_registry import ureg
import numpy as np
import matplotlib.pyplot as plt


def simple_setup(input_parameter, height,fc,load):

    out = check_beam_design(span=input_parameter['beamExSpan']*ureg(input_parameter['beamExSpanUnit']),
                            width=input_parameter['beamExWidth']*ureg(input_parameter['beamExWidthUnit']),
                            height=height,
                            point_load=load,
                            distributed_load=0*ureg('N/mm'),
                            compr_str_concrete=fc,
                            yield_str_steel=input_parameter['beamExYieldStrSteel']*ureg(input_parameter['beamExYieldStrSteelUnit']),
                            steel_dia_bu=input_parameter['beamExSteelDiaBu']*ureg(input_parameter['beamExSteelDiaBuUnit']),
                            cover_min=input_parameter['beamExCoverMin']*ureg(input_parameter['beamExCoverMinUnit']))

    return out


def beam_design_plot(input_parameter, n, fig_path:str = "beam_design_plot.pdf"):
    # DISCLAMER ;)
    # this is a wip, to see if it works as expected.
    # this will be made pretty etc. in my paper branch, once this is merged
    # TODO:
    #       integrate this into dodo and tex

    min_fc = 1
    max_fc = 100
    fc_unit = 'N/mm^2'
    min_load = 10
    max_load = 200
    load_unit = 'kN'
    min_height = 150
    max_height = 500
    height_unit = 'mm'
    crosssection_unit = 'cm^2'

    def get_list(min,max,n):
        step = (max-min)/(n-1)
        return np.arange(min, max+step, step)

    def get_plot_lists(order_list, values_dict, input_parameter):
        print(f'compute: {order_list[0]} - {order_list[1]}')
        x_list = values_dict[order_list[0]]['list']
        y_list = values_dict[order_list[1]]['list']
        constant = values_dict[order_list[2]]['constant']


        order_list = [0 if item == 'height' else item for item in order_list]
        order_list = [1 if item == 'fc' else item for item in order_list]
        order_list = [2 if item == 'load' else item for item in order_list]

        crosssections = np.zeros(shape=(len(x_list), len(y_list)))
        fc_errors = np.zeros(shape=(len(x_list), len(y_list)))
        A_errors = np.zeros(shape=(len(x_list), len(y_list)))

        values = [0,0,0]
        values[order_list[2]] = constant
        for i,x in enumerate(x_list):
            values[order_list[0]] = x
            for j,y in enumerate(y_list):
                values[order_list[1]] = y

                out = simple_setup(input_parameter,*values)
                crosssection = out['crosssection']
                crosssection.ito(crosssection_unit)
                crosssection = crosssection.magnitude
                crosssections[i][j] = crosssection
                fc_errors[i][j] = out['fc_error']
                A_errors[i][j] = out['A_error']

        return crosssections, fc_errors, A_errors

    def plot_contour(ax,x,y,Z,xlabel,ylabel,title,colorbar=False,max=None):
        if not max:
            max = np.nanmax(Z)

        min = np.nanmin(Z)
        # # TODO: testing...
        # max = 16.08

        n_levels = 8
        level_step = max/(n_levels-1)
        levels = np.arange(min,max+level_step,level_step)

        # # TODO: testing:
        # levels = np.append(levels, [200])

        X, Y = np.meshgrid(y, x)
        #cs = ax.contourf(X, Y, Z, vmin=0, vmax=max-level_step,levels=levels)
        cs = ax.contourf(X, Y, Z, vmin=min, vmax=max-level_step, levels=levels)
        ax.set_ylabel(xlabel)
        ax.set_xlabel(ylabel)
        ax.set_title(title)

        if colorbar:
            cbar = plt.colorbar(cs)
            cbar.set_label(f'total steel cross section in {ureg(crosssection_unit).units}')


    constant_fc = input_parameter['beamExComprStrConcreteC'] * ureg(input_parameter['beamExComprStrConcreteCUnit'])
    constant_height = input_parameter['beamExHeightC'] * ureg(input_parameter['beamExHeightCUnit'])
    constant_load = input_parameter['beamExPointLoadC'] * ureg(input_parameter['beamExPointLoadCUnit'])

    # generate lists
    fc_list = get_list(min_fc, max_fc, n) * ureg(fc_unit)
    load_list = get_list(min_load,max_load,n) * ureg(load_unit)
    height_list = get_list(min_height,max_height,n) * ureg(height_unit)

    values = {'fc' : {'list' : get_list(min_fc, max_fc, n) * ureg(fc_unit),
                      'constant' : input_parameter['beamExComprStrConcreteC'] * ureg(input_parameter['beamExComprStrConcreteCUnit'])},
              'load': {'list': get_list(min_load,max_load,n) * ureg(load_unit),
                       'constant': input_parameter['beamExPointLoadC'] * ureg(input_parameter['beamExPointLoadCUnit'])},
              'height': {'list': get_list(min_height,max_height,n) * ureg(height_unit),
                         'constant': input_parameter['beamExHeightC'] * ureg(input_parameter['beamExHeightCUnit'])}}

    chosen_plots = [('fc', 'load', 'height'),
                    ('height', 'load', 'fc'),
                    ('height', 'fc', 'load')]

    # empty lists
    crosssection_plots = [0]*len(chosen_plots)
    fc_error_plots = [0]*len(chosen_plots)
    A_error_plots = [0]*len(chosen_plots)
    plot_values = {}
    for i, plot in enumerate(chosen_plots):
        crosssection_plots[i], fc_error_plots[i], A_error_plots[i] = get_plot_lists(plot, values, input_parameter)
        plot_values[plot] = {'crosssections' : crosssection_plots[i],
                             'fc_errors' : fc_error_plots[i],
                             'A_errors' : A_error_plots[i]}




    fig, axs = plt.subplots(3, 3, figsize=(15, 12))


    #max_crossection = max([np.nanmax(crossections_1),np.nanmax(crossections_2),np.nanmax(crossections_3)])
    #    # plot figure 1

    for i, plot in enumerate(chosen_plots):
        plot_contour(axs[0][i], values[plot[0]]['list'], values[plot[1]]['list'], plot_values[plot]['crosssections'],
                     xlabel=f'{plot[0]} in {values[plot[0]]["list"].units}',
                     ylabel=f'{plot[1]} in {values[plot[1]]["list"].units}',
                     title=f'{plot[2]} = {values[plot[2]]["constant"]}', colorbar=True)

        plot_contour(axs[1][i], values[plot[0]]['list'], values[plot[1]]['list'], plot_values[plot]['fc_errors'],
                     xlabel=f'{plot[0]} in {values[plot[0]]["list"].units}',
                     ylabel=f'{plot[1]} in {values[plot[1]]["list"].units}',
                     title=f'{plot[2]} = {values[plot[2]]["constant"]}', colorbar=True)

        plot_contour(axs[2][i], values[plot[0]]['list'], values[plot[1]]['list'], plot_values[plot]['A_errors'],
                     xlabel=f'{plot[0]} in {values[plot[0]]["list"].units}',
                     ylabel=f'{plot[1]} in {values[plot[1]]["list"].units}',
                     title=f'{plot[2]} = {values[plot[2]]["constant"]}', colorbar=True)

    fig.tight_layout()
    fig.savefig(fig_path)




if __name__ == "__main__":

    beam_design_example_parameters = {
        'beamExSpan': 675,
        'beamExSpanUnit': 'cm',
        'beamExWidth': 200,
        'beamExWidthUnit': 'mm',
        'beamExYieldStrSteel': 500,
        'beamExYieldStrSteelUnit': 'N/mm^2',
        'beamExSteelDiaBu': 10,
        'beamExSteelDiaBuUnit': 'mm',
        'beamExCoverMin': 2.5,
        'beamExCoverMinUnit': 'cm',
        'beamExHeightC': 450,
        'beamExHeightCUnit': 'mm',
        'beamExPointLoadC': 50,
        'beamExPointLoadCUnit': 'kN',
        'beamExComprStrConcreteC': 40,
        'beamExComprStrConcreteCUnit': 'N/mm^2',
    }

    beam_design_plot(beam_design_example_parameters, n=10)