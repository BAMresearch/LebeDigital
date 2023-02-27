from lebedigital.unit_registry import ureg
from lebedigital.demonstrator_scripts.kpi_from_fem import kpi_from_fem
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx
import pytest
import pint_pandas
import pint
import pandas as pd


def test_kpi_from_fem():

    pint.set_application_registry(ureg)  # required to use the same registry
    df = pd.DataFrame({
        "time": pd.Series([0,10,20], dtype="pint[hours]"),
        "temperature": pd.Series([10,40,80], dtype="pint[degree_Celsius]"),
        "yield": pd.Series([-40,-20,20], dtype="pint[dimensionless]"),
    })

    Q_ = ureg.Quantity
    limit = Q_(70,ureg.degC)

    results = kpi_from_fem(df, limit)

    assert results['max_reached_temperature'].magnitude == 80
    assert results['time_max_reached_temperature'].magnitude == 20
    assert results['check_reached_temperature'].magnitude == -10
    assert results['time_of_demolding'].magnitude == 15
