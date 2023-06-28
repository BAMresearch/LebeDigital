import numpy as np
import pandas as pd
import pint
import pint_pandas
from scipy.optimize import curve_fit

from lebedigital.unit_registry import ureg


def kpi_from_fem(df, limit_temp):
    """
    compute KPIs from simulation output

    the maximum reached temperature and the time is extracted and the difference to the limit temperature given
    the time when the yield value equals zero is approximated
        - a yield surface is defined as the elastic limit, based on von Mises and Rankine criteria
        - for each timestep the maximum value is given (positive number indicate a exceeding of the elastic range)
        - once a timestep returns a negative value as maximum, the beam is considered to not fail
        - from this the time can be linearly interpolated then the max value is approximately zero

    Parameters
    ----------
    df : panda-pint dataframe
        data frame with the columns "time","temperature","yield"
    limit_temp: float in pint units
        maximum allowed temperature

    Returns
    -------
    results : dict with pint unit
        - 'max_reached_temperature' in degree_Celsius
        - 'time_max_reached_temperature' in hours
        - 'check_reached_temperature' in degree_Celsius
        - 'time_of_demolding' in hours
    """
    # initialze dictionary
    results = {}

    # Part 1 : interpolate time when yield == 0
    # this is based on the assumption that the sign of the max yield value only changes once.
    # if there are cases where that is not ture, this approach need to be improved.
    # TODO: catch multiple possible zero values, minimum: throw an error if this happens

    # make sure the columns have expected units
    df["time"] = df["time"].pint.to("seconds")
    df["temperature"] = df["temperature"].pint.to("degree_Celsius")
    df["yield"] = df["yield"].pint.to("dimensionless")

    # subsequent operations have problems with the units (interpolation) therefore convert to standard df
    df = df.pint.dequantify()

    # adding new row with zero yield
    new_line = pd.DataFrame(
        {("time", "second"): [np.nan], ("temperature", "degree_Celsius"): [np.nan], ("yield", "dimensionless"): [0.0]}
    )
    df = pd.concat([df, new_line], ignore_index=True)

    # sort table
    df = df.sort_values(by=[("yield", "dimensionless")], ascending=False)

    # interpolate missing values
    df = df.interpolate(method="linear", limit_direction="forward")

    # locate time of demoldung
    results["time_of_demolding"] = df.at[
        df.loc[df[("yield", "dimensionless")] == 0.0].index.values[0], ("time", "second")
    ] * ureg("s")

    # check if the interpolation worked
    if np.isnan(results["time_of_demolding"]):
        # based on https://stackoverflow.com/questions/22491628/extrapolate-values-in-pandas-dataframe
        # extrapolate missing values

        # Function to curve fit to the data
        def func(x, a, b, c):
            return a * (x**2) + b * x + c

        # Initial parameter guess, just to kick off the optimization
        guess = (0.5, 0.5, 0.5)

        # Create copy of data to remove NaNs for curve fitting
        fit_df = df.dropna()

        # Place to store function parameters for each column
        col_params = {}

        # Curve fit each column
        for col in fit_df.columns:
            # Get x & y
            x = fit_df.index.astype(float).values
            y = fit_df[col].values
            # Curve fit column and get curve parameters
            params = curve_fit(func, x, y, guess)
            # Store optimized parameters
            col_params[col] = params[0]

        # Extrapolate each column
        for col in df.columns:
            # Get the index values for NaNs in the column
            x = df[pd.isnull(df[col])].index.astype(float).values
            # Extrapolate those points with the fitted function
            df[col][x] = func(x, *col_params[col])

        results["time_of_demolding"] = df.at[
            df.loc[df[("yield", "dimensionless")] == 0.0].index.values[0], ("time", "second")
        ] * ureg("s")

    # changing units, because we can
    results["time_of_demolding"].ito("h")

    # requantify the dataframe
    df = df.pint.quantify()

    # Part 2 : get max temp, with time -> compute difference to limit
    results["max_reached_temperature"] = df["temperature"].max()

    results["time_max_reached_temperature"] = df.loc[df["temperature"] == results["max_reached_temperature"]][
        "time"
    ].tolist()[0]

    # changing units, because we can
    results["time_max_reached_temperature"].ito("hours")

    # difference between limit temp and reached maximum
    results["check_reached_temperature"] = limit_temp - results["max_reached_temperature"]

    return results
