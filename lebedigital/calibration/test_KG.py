import rdflib
import pandas as pd
import numpy as np
import fenics_concrete

def cylinder_simulation(parameters):
        """running the cylinder simulation
        Optimized for a linear problem as we only need to compute once

        Parameter
        ---------
        parameters: dic
                dictionary with simulation parameters

        Return
        ------
        slope: float
                slope of the linear load displacement function
                """

        parameters['nu'] = 0.25
        parameters['dim'] = 3

        parameters['mesh_density'] = 10
        parameters['log_level'] = 'WARNING'
        parameters['bc_setting'] = 'fixed'

        test_load = -0.05

        experiment = fenics_concrete.ConcreteCylinderExperiment(parameters)

        problem = fenics_concrete.LinearElasticity(experiment, parameters)
        sensor = fenics_concrete.sensors.ReactionForceSensorBottom()

        problem.add_sensor(sensor)
        problem.experiment.apply_displ_load(test_load)
        problem.solve()  # solving this
        measured_force = problem.sensors[sensor.name].data[-1]

        # compute slope of linear problem
        slope = measured_force/test_load

        return slope,measured_force

def query_objects(queries, graph):
        # function to get objects from specific subject, predicate pairs, knowing there is only one result
        results = {}
        for query in queries:
                q = f"""
                        SELECT ?object
                        WHERE {{
                                {queries[query]['subject']}  {queries[query]['predicate']} ?object
                        }}
                """

                result = graph.query(q)
                for r in result:
                        results[query] = r["object"]

        return results

# define queries
queries = {
        'diameter' :  {'subject' : 'ns1:informationbearingentity1', 'predicate' : 'ns1:has_decimal_value' },
        'length' :    {'subject' : 'ns1:informationbearingentity2', 'predicate' : 'ns1:has_decimal_value' },
        'path' :      {'subject' : 'ns1:informationbearingentity9', 'predicate' : 'ns9:has_text_value' },
        'file_name' : {'subject' : 'ns1:informationbearingentity8', 'predicate' : 'ns9:has_text_value' }
}

# Path to KG, for testing
path_to_KG = '../../usecases/MinimumWorkingExample/emodul/knowledge_graphs/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.tll'

# initialize the graph
knowledge_graph = rdflib.Graph()
knowledge_graph.parse(path_to_KG, format='ttl')

# queries
results = query_objects(queries,knowledge_graph)

print('STEP 1 : getting data from knowledge graph')
# geometry
length = float(results['length'])
print('Length: ', length)
diameter = float(results['diameter'])
print('Diameter: ',diameter)

# experimental data
file_path = results['path'] + '/' + results['file_name']
print('File path: ',file_path)

print('STEP 2 : reading processed data')
# read csv into lists
df = pd.read_csv(file_path)

# trying to get rid of stupid values, i.e. positive displacements and forces
# TODO: is there anything else ot "clean" the raw data?
df = df[df['Transducer 1[mm]'] < 0]
df = df[df['Transducer 2[mm]'] < 0]
df = df[df['Transducer 3[mm]'] < 0]
df = df[df['Force [kN]'] < 0]

force_list = df['Force [kN]'].to_list()
displacement_list_1 = np.array(df['Transducer 1[mm]'].to_list())
displacement_list_2 = np.array(df['Transducer 2[mm]'].to_list())
displacement_list_3 = np.array(df['Transducer 3[mm]'].to_list())

print('STEP 3 : run a simulation')
# testing the cylinder simulation
parameters = fenics_concrete.Parameters()  # using the current default values

parameters['radius'] = diameter/2  # in mm
parameters['height'] = length   # in mm
parameters['E'] = 100  # in kN/mm^2
print(' - Youngs modulus: ', parameters['E'])

simulation_slope, force = cylinder_simulation(parameters)

print('Slope:', simulation_slope)

force_list_1 = displacement_list_1*simulation_slope
force_list_2 = displacement_list_2*simulation_slope
force_list_3 = displacement_list_3*simulation_slope

print('STEP 4 : TODO optimization')


#
#



