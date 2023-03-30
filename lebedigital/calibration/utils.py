## -- TUM PKM --- atul.agrawal@tum.de--- ##
import matplotlib.pyplot as plt
import rdflib
import pandas as pd
import numpy as np
import os
from pathlib import Path
import sys
import fenics_concrete

baseDir1 = Path(__file__).resolve().parents[1]
sys.path.append(os.path.join(os.path.join(baseDir1, "knowledgeGraph"), "emodul"))
sys.path.append(os.path.join(baseDir1, "Data"))


# import emodul_query

def read_exp_data_E_mod(path: str, exp_name: str,length:float,diameter:float) -> dict:
    """
    Reads in the experiment data for a specified experiment and stores the extracted results in the dict.
    The arguments to the provided by the Knowledge Graph.

    Parameters
    ----------
    path : str Path to the folder where the experimental data in .csv format is stored.
    exp_name : str the experiment name.csv
    length : float the length of the specimen
    diameter : the diameter of the specimen

    Returns
    -------
    output : dict
    values with the keys : exp_name, height, diameter, displacement, force
    """
    file_path = os.path.join(path, exp_name)
    results = {}
    results["length"] = length
    results["diameter"] = diameter
    df = extract_third_load_cycle(data_path=file_path)
    df["displacement"] = (
        df["Transducer 1[mm]"] + df["Transducer 2[mm]"] + df["Transducer 3[mm]"]
    ) / 3
    # df['stress'] = df['Force [kN]'] / (np.pi * (float(results['diameter']) / 2) ** 2)

    output = {
        "exp_name": exp_name,
        "height": float(results["length"]),
        "diameter": float(results["diameter"]),
        "displacement": np.array(df["displacement"]),
        "force": np.array(df["Force [kN]"]),
    }
    return output


class PosteriorPredictive:
    def __init__(
        self, forward_solver, known_input_solver, parameter=None
    ):

        """

        Parameters
        ----------
        forward_solver :
            Should be wrapped in a function which takes in two args. arg1 = random input, arg2 = known input
        parameter :
            An array of samples. Currently parameters dont support a distrutbion, just MCMC samples. To get an approx.
            for parameter dist, we need VI based posterior approx
        known_input_solver :
            Pass known input value to the solver in whichever format the forward_solve object needs.
        """

        # assert (parameter == np.array), "Parameter needs to be np.ndarray"
        # assert parameter.shape[1] == 1, "Currently support only 1 parameter"

        self._forward_solver = forward_solver
        if parameter is not None:
            self._parameter = parameter
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
        raise NotImplementedError("...")


def extract_third_load_cycle(
    data_path: str, threshold=1, vizualize=False
) -> pd.DataFrame:
    """
    Script to extract third loading cycle of the load-displacement curve for a given experiment data as a .csv file

    Parameters
    ----------
    data_path : The path to the experimental data stored as a .csv file.
    threshold : (not recommended to be modified)
    vizualize : To viz. the original data and the extracted data
    Returns
    -------
    data_third_loading : Dataframe containing the third loading cycle
    """
    # Load data from .csv file
    data = pd.read_csv(data_path, skipfooter=5)

    # Extract indices where there is a change in slope
    slope_2 = np.diff(
        data["Force [kN]"], n=2
    )  # double diff to identify the sharp points
    change_indices = np.where(np.abs(slope_2) > threshold)[0] + 2  # Finding the indices

    ## drop the indices which are close together.
    change_indices_filtered = []
    tolerance = 8  # In hopes this will cover all the edge cases
    for i, value in enumerate(change_indices):
        if i == 0 or abs(value - change_indices[i - 1]) > tolerance:
            change_indices_filtered.append(value)

    # Select indices for the third loading cycle and update the data
    idx = [
        change_indices_filtered[-4],
        change_indices_filtered[-3],
    ]  # skipping the last two change in slopes
    data_third_loading = data.loc[idx[0] : idx[1]]

    # Plot and see the data
    if vizualize:
        # plot original data
        plt.plot(data["Force [kN]"], "r", label="Original Data")
        plt.plot(data_third_loading["Force [kN]"], "g", label="Third load cycle")
        plt.legend()
        plt.show()
    return data_third_loading
