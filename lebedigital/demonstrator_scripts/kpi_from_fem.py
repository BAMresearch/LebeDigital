from lebedigital.unit_registry import ureg
import numpy as np
import pint_pandas
import pint
import pandas as pd


def kpi_from_fem(df,limit_temp):
    """
    compute KPIs from simulation output

    the maximum reached temperature and the time is extracted and the difference to the limit temperature given
    the time when the yield value equals zero is approximated

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

    # Part 1 : get max temp, with time -> compute difference to limit

    results['max_reached_temperature'] = df['temperature'].max()
    results['time_max_reached_temperature'] = \
        df.loc[df['temperature'] == results['max_reached_temperature']]['time'].tolist()[0]

    # changing units, because we can
    results['time_max_reached_temperature'].ito('hours')

    # difference between limit temp and reached maximum
    results['check_reached_temperature'] = limit_temp - results['max_reached_temperature']

    # Part 2 : interpolate time when yield == 0

    # make sure the columns have expected units
    df['time'] = df['time'].pint.to("seconds")
    df['temperature'] = df['temperature'].pint.to("degree_Celsius")
    df['yield'] = df['yield'].pint.to("dimensionless")

    # subsequent operations have problems with the units (interpolation) therefore convert to standard df
    df = df.pint.dequantify()

    # adding new row with zero yield
    new_line = pd.DataFrame({('time', 'second'): [np.nan], ('temperature', 'degree_Celsius'): [np.nan], ('yield', 'dimensionless'): [0.0]})
    df = pd.concat([df, new_line],  ignore_index=True)

    # sort table
    df = df.sort_values(by=[('yield', 'dimensionless')], ascending=False)

    # interpolate missing values
    df = df.interpolate(method='linear', limit_direction='forward')

    # locate time of demoldung
    results['time_of_demolding'] = df.at[df.loc[df[('yield', 'dimensionless')] == 0.0].index.values[0],('time', 'second')] * ureg('s')

    # changing units, because we can
    results['time_of_demolding'].ito('h')

    return results


if __name__ == "__main__":
    # test while developing this

    pint.set_application_registry(ureg)  # required to use the same registry
    df = pd.DataFrame({
        "time": pd.Series([0,10,20], dtype="pint[hours]"),
        "temperature": pd.Series([10,40,80], dtype="pint[degree_Celsius]"),
        "yield": pd.Series([-40,-20,20], dtype="pint[dimensionless]"),
    })

    Q_ = ureg.Quantity
    limit = Q_(70,ureg.degC)

    print(kpi_from_fem(df, limit))
