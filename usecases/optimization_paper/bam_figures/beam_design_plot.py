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

    constant_fc = input_parameter['beamExComprStrConcreteC'] * ureg(input_parameter['beamExComprStrConcreteCUnit'])
    constant_height = input_parameter['beamExHeightC'] * ureg(input_parameter['beamExHeightCUnit'])
    constant_load = input_parameter['beamExPointLoadC'] * ureg(input_parameter['beamExPointLoadCUnit'])

    min_fc = 1
    max_fc = 100
    fc_unit = 'N/mm^2'
    min_load = 10
    max_load = 200
    load_unit = 'kN'
    min_height = 150
    max_height = 500
    height_unit = 'mm'





    def get_list(min,max,n):
        step = (max-min)/(n-1)
        return np.arange(min, max+step, step)

    fc_list = get_list(min_fc, max_fc, n) * ureg(fc_unit)
    load_list = get_list(min_load,max_load,n) * ureg(load_unit)
    height_list = get_list(min_height,max_height,n) * ureg(height_unit)

    if True:
        print('compute: fc - load')
        crossections_1 = np.zeros(shape=(len(fc_list), len(load_list)))

        for i,fc in enumerate(fc_list):
            for j,load in enumerate(load_list):
                try:
                    out = simple_setup(input_parameter, constant_height,fc ,load )
                    crosssection = out['crosssection']
                    crosssection.ito('cm^2')
                    crosssection = crosssection.magnitude
                except:
                    crosssection = np.nan
                crossections_1[i][j] = crosssection

    if True:
        print('compute: height - load')
        # plot height vs load
        crossections_2 = np.zeros(shape=(len(height_list), len(load_list)))

        for i,height in enumerate(height_list):
            for j,load in enumerate(load_list):
                try:
                    out = simple_setup(input_parameter, height,constant_fc,load )
                    crosssection = out['crosssection']
                    crosssection.ito('cm^2')
                    crosssection = crosssection.magnitude
                except:
                    crosssection = np.nan
                crossections_2[i][j] = crosssection


    if True:
        print('compute: height - fc')
        #plot fc vs height

        crossections_3 = np.zeros(shape=(len(height_list), len(fc_list)))

        for i,height in enumerate(height_list):
            for j,fc in enumerate(fc_list):
                try:
                    out = simple_setup(input_parameter, height,fc,constant_load)
                    crosssection = out['crosssection']
                    crosssection.ito('cm^2')
                    crosssection = crosssection.magnitude
                except:
                    crosssection = np.nan
                crossections_3[i][j] = crosssection



        def plot_contour(ax,x,y,Z,xlabel,ylabel,title,colorbar=False,max=None):
            if not max:
                max = np.nanmax(Z)

            n_levels = 8
            level_step = max/(n_levels-1)
            levels = np.arange(0,max+level_step,level_step)

            X, Y = np.meshgrid(y, x)
            cs = ax.contourf(X, Y, Z, vmin=0, vmax=max-level_step,levels=levels)
            ax.set_ylabel(xlabel)
            ax.set_xlabel(ylabel)
            ax.set_title(title)



            if colorbar:
                cbar = plt.colorbar(cs)
                cbar.set_label('cross section')


        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        max_crossection = max([np.nanmax(crossections_1),np.nanmax(crossections_2),np.nanmax(crossections_3)])
        # plot figure 1
        plot_contour(axs[2],fc_list,load_list,crossections_1,
                     xlabel=f'fc in {fc_list.units}',
                     ylabel=f'load in {load_list.units}',
                     title=f'height = {constant_height}',max=max_crossection, colorbar=True)
        plot_contour(axs[1],height_list,load_list,crossections_2,
                     xlabel=f'height in {height_list.units}',
                     ylabel=f'load in {load_list.units}',
                     title=f'fc = {constant_fc}',max=max_crossection)
        plot_contour(axs[0],height_list,fc_list,crossections_3,
                     xlabel=f'height in {height_list.units}',
                     ylabel=f'fc in {fc_list.units}',
                     title=f'load = {constant_load}',max=max_crossection)

        # plt.subplots_adjust(wspace=0.4)
        fig.tight_layout()

        fig.savefig(fig_path)
        #plt.show()




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