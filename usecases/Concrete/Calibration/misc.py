## -- TUM PKM --- atul.agrawal@tum.de--- ##

import pandas as pd
import numpy as np
from usecases.Concrete.knowledgeGraph.query import input_emodul_data_for_calibration


def load_experimental_data(exp_name, skip_init, skip_last):
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

    df = pd.read_csv(input_emodul_data_for_calibration(exp_name)['processedDataPath'], skipfooter=skip_last)
    df = df.drop(labels=range(0, skip_init), axis=0)

    df['displacement'] = (df['Transducer 1[mm]'] + df['Transducer 2[mm]'] + df['Transducer 3[mm]']) / 3
    dia = input_emodul_data_for_calibration(exp_name)['specimenDiameter']
    df['stress'] = df['Force [kN]'] / (np.pi * (dia / 2) ** 2)
    height = input_emodul_data_for_calibration(exp_name)['specimenLength']

    output = {
        'height': height,
        'diameter': dia,
        'displacement': np.array(df['displacement']),
        'stress': np.array(df['stress'])}
    return output
