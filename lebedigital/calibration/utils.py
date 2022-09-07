## -- TUM PKM --- atul.agrawal@tum.de--- ##

import rdflib
import pandas as pd
import numpy as np
import os
from pathlib import Path
import sys

baseDir1 = Path(__file__).resolve().parents[1]
sys.path.append(os.path.join(os.path.join(baseDir1, 'knowledgeGraph'), 'emodul'))
sys.path.append(os.path.join(baseDir1, 'Data'))

#import emodul_query

def query_KG(path: str, skip_last = 145, skip_init = 330 ) -> dict:
    """
    Queries KG, the updated one done by Illias (I think from KG team). Lot of funky words as I dont know what they are.
    Returns:
    Args:
        path ():
        skip_last (): Last datavalues to skip as we only need the third loading cycle.
        skip_init (): Initial datavalues to skip as we only need the third loading cycle.

    Returns:

    """
    """
    Queries KG, the updated one done by Illias (I think from KG team). Lot of funky words as I dont know what they are.
    Returns:

    """

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
        'diameter': {'subject': 'ns1:informationbearingentity1', 'predicate': 'ns1:has_decimal_value'},
        'length': {'subject': 'ns1:informationbearingentity2', 'predicate': 'ns1:has_decimal_value'},
        'path': {'subject': 'ns1:informationbearingentity9', 'predicate': 'ns9:has_text_value'},
        'file_name': {'subject': 'ns1:informationbearingentity8', 'predicate': 'ns9:has_text_value'}
    }

    # Path to KG, for testing
    #path_to_KG = '../../usecases/MinimumWorkingExample/emodul/knowledge_graphs/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4.ttl'
    path_to_KG = path
    # initialize the graph
    knowledge_graph = rdflib.Graph()
    knowledge_graph.parse(path_to_KG, format='ttl')

    # queries
    results = query_objects(queries, knowledge_graph)

    # read csv into dataframe
    file_path = results['path'] + '/' + results['file_name']
    #skip_last = 145 # hard coded to get the third loading cycle. need somthing better
    #skip_init = 330
    df = pd.read_csv(file_path, skipfooter=skip_last, engine='python')
    df = df.drop(labels=range(0, skip_init), axis=0)
    df['displacement'] = (df['Transducer 1[mm]'] + df['Transducer 2[mm]'] + df['Transducer 3[mm]']) / 3
    df['stress'] = df['Force [kN]'] / (np.pi * (float(results['diameter']) / 2) ** 2)

    output = {
        'height': float(results['length']),
        'diameter': float(results['diameter']),
        'displacement': np.array(df['displacement']),
        'stress': np.array(df['stress'])}
    return output

def load_experimental_data(exp_name, skip_init, skip_last, KG=False, path=None):
    """
    Function to interact with the knowledge graph database and get experimental data
    # TODO : Clarify the units used. Now Force :kN, dispalcement : mm
    Args:
        skip_init (int): Initial datavalues to skip as we only need the third loading cycle.
        skip_last (int): Last datavalues to skip as we only need the third loading cycle.
        exp_name (string): The experiment name

    Returns:
        output (dict):
            The experiment values with keys as 'height' (int), 'diameter' (int), 'displacement' (np.array), 'stress' (np.array)
    """
    if KG is True:
        df = pd.read_csv(emodul_query.input_emodul_data_for_calibration(exp_name)['processedDataPath'],
                         skipfooter=skip_last)
    if KG is not True:  # manually passing the processed data path
        df = pd.read_csv(path, skipfooter=skip_last)
    df = df.drop(labels=range(0, skip_init), axis=0)

    if KG is True:
        dia = emodul_query.input_emodul_data_for_calibration(exp_name)['specimenDiameter']
        height = emodul_query.input_emodul_data_for_calibration(exp_name)['specimenLength']
    if KG is not True:
        dia = 98.6
        height = 300.2

    df['displacement'] = (df['Transducer 1[mm]'] + df['Transducer 2[mm]'] + df['Transducer 3[mm]']) / 3
    df['stress'] = df['Force [kN]'] / (np.pi * (dia / 2) ** 2)

    output = {
        'height': height,
        'diameter': dia,
        'displacement': np.array(df['displacement']),
        'stress': np.array(df['stress'])}
    return output

class PosteriorPredictive:
    """
    Will be depricated when it is integrated in probeye
    """

    def __init__(self, forward_solver, known_input_solver, parameter = None, query_kg = None):

        """
        TODO: Add option to add point estimates of parameters if forward solver is expensive

        Parameters
        ----------
        forward_solver :
            Defines the forward model. Check out forward_model.py to see a template for
            the forward model definition if needed. The user will then have to derive his own
            forward model from that base class. Examples can be found in the package
            directory tests/integration_tests. user can also specify his/her own forward model,
            which just takes in input.
            Should be wrapped in a fucntion which takes in two args. arg1 = random input, arg2 = known input
        parameter :
            An array of samples. Currently parameters dont support a distrutbion, just MCMC samples. To get an approx.
            for parameter dist, we need VI based posterior approx
        known_input_solver :
            Pass known input value to the solver in whichever format the forward_solve object needs.
        query_kg :
            A query script (with appropriate arguments) which outputs the sample or the stats of the inferred parameter
            needed for the forward model
        """

        #assert (parameter == np.array), "Parameter needs to be np.ndarray"
        #assert parameter.shape[1] == 1, "Currently support only 1 parameter"

        self._forward_solver = forward_solver
        if parameter is not None:
            self._parameter = parameter
        if query_kg is not None:
            # run the query script
            para = query_kg()
            self._parameter = para
        self._known_input = known_input_solver
        self._mean = None
        self._std = None
        self._samples = None

    def get_stats(self, samples: int) -> tuple:
        """
        Returns mean and s.d of the posterior predictive. Simple Monte Carlo based approximation

        Parameters
        ----------
        samples :

        Returns
        -------

        """
        # Monte carlo step
        output = []
        for i in range(0, samples):
            y = self._forward_solver(self._parameter[i], self._known_input)
            output.append(y)

        # get the posterior pred stats
        mean = np.mean(output, axis=0)
        sd = np.std(output, axis=0)

        self._mean = mean
        self._std = sd
        self._samples = output
        return mean, sd

    def plot_posterior_predictive(self):
        raise NotImplementedError(
            "...")
