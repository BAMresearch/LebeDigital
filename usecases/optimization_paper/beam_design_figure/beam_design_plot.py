from lebedigital.demonstrator_scripts.beam_design import check_beam_design
from lebedigital.unit_registry import ureg
import numpy as np
import matplotlib.pyplot as plt


def simple_setup(height,fc,load):

    out = check_beam_design(span=6750*ureg('mm'),
                            width=200*ureg('mm'),
                            height=height*ureg('mm'),
                            point_load=load*ureg('kN'),
                            distributed_load=0*ureg('N/mm'),
                            compr_str_concrete=fc*ureg('N/mm^2'),
                            yield_str_steel=500*ureg('N/mm^2'),
                            steel_dia_bu=12*ureg('mm'),
                            cover_min=2.5*ureg('cm'))

    return out


def beam_design_plot():
    # DISCLAMER ;)
    # this is a wip, to see if it works as expected.
    # this will be made pretty etc. in my paper branch, once this is merged
    # TODO: add units to plots
    #       add input to a file
    #       integrate this into dodo and tex

    def get_list(min,max,n):
        step = (max-min)/(n-1)
        return np.arange(min, max+step, step)

    n = 10

    if True:
        print('compute: fc - load')
        fc_list_1 = get_list(1,100,n)
        load_list_1 = get_list(10,200,n)
        crossections_1 = np.zeros(shape=(len(fc_list_1), len(load_list_1)))
        height = 450

        for i,fc in enumerate(fc_list_1):
            for j,load in enumerate(load_list_1):
                try:
                    out = simple_setup(height,fc,load)
                    crosssection = out['crosssection']
                    crosssection.ito('cm^2')
                    crosssection = crosssection.magnitude
                except:
                    crosssection = np.nan
                crossections_1[i][j] = crosssection


    if True:
        print('compute: height - load')
        # plot height vs load
        height_list_2 = get_list(130,480,n)
        load_list_2 = get_list(10,200,n)
        crossections_2 = np.zeros(shape=(len(height_list_2), len(load_list_2)))
        fc = 50

        for i,height in enumerate(height_list_2):
            for j,load in enumerate(load_list_2):
                try:
                    out = simple_setup(height,fc,load)
                    crosssection = out['crosssection']
                    crosssection.ito('cm^2')
                    crosssection = crosssection.magnitude
                except:
                    crosssection = np.nan
                crossections_2[i][j] = crosssection


    if True:
        print('compute: height - fc')
        #plot fc vs height
        height_list_3 = get_list(170,790,n)
        fc_list_3 = get_list(1,100,n)
        load = 40
        crossections_3 = np.zeros(shape=(len(height_list_3), len(fc_list_3)))

        for i,height in enumerate(height_list_3):
            for j,fc in enumerate(fc_list_3):
                try:
                    out = simple_setup(height,fc,load)
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
        plot_contour(axs[2],fc_list_1,load_list_1,crossections_1,
                     xlabel='fc',ylabel='load',title='height=450',max=max_crossection, colorbar=True)
        plot_contour(axs[1],height_list_2,load_list_2,crossections_2,
                     xlabel='height',ylabel='load',title='fc=50',max=max_crossection)
        plot_contour(axs[0],height_list_3,fc_list_3,crossections_3,
                     xlabel='height',ylabel='fc',title='load=40',max=max_crossection)
        plt.show()






if __name__ == "__main__":
    # homogenization_example_parameters = {
    #     'pasteE': 30e9,
    #     'pasteEunit': 'Pa',
    #     'pastenu': 0.2,
    #     'pastenuunit': 'dimensionless',
    #     'pasteC': 870,
    #     'pasteCunit': 'J/kg/K',
    #     'pastekappa': 1.8,
    #     'pastekappaunit': 'W/m/K',
    #     'pasterho': 2400,
    #     'pasterhounit': 'kg/m^3',
    #     'pastefc': 30e6,
    #     'pastefc_unit': 'Pa',
    #     'pasteQ': 250000,
    #     'pasteQunit': 'J/kg',
    #     'aggregatesE': 25e9,
    #     'aggregatesE_unit': 'Pa',
    #     'aggregatesnu': 0.3,
    #     'aggregatesnu_unit': 'dimensionless',
    #     'aggregatesC': 840,
    #     'aggregatesCunit': 'J/kg/K',
    #     'aggregateskappa': 0.8,
    #     'aggregateskappaunit': 'W/m/K',
    #     'aggregatesrho': 2600,
    #     'aggregatesrhounit': 'kg/m^3'
    # }
    beam_design_plot()

    #create_homogenization_figure(homogenization_example_parameters)