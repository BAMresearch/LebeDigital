import pandas as pd
import pint
import pint_pandas
import pytest
from pint.testsuite.helpers import assert_quantity_almost_equal as assert_approx

from lebedigital.demonstrator_scripts.kpi_from_fem import kpi_from_fem
from lebedigital.unit_registry import ureg


#
def test_kpi_from_fem():
    pint.set_application_registry(ureg)  # required to use the same registry
    df = pd.DataFrame(
        {
            "time": pd.Series([0, 10, 20], dtype="pint[hours]"),
            "temperature": pd.Series([10, 40, 80], dtype="pint[degree_Celsius]"),
            "yield": pd.Series([40, 20, -20], dtype="pint[dimensionless]"),
        }
    )

    Q_ = ureg.Quantity
    limit_T = Q_(70, ureg.degC)
    limit_time = 6 * ureg.hours

    results = kpi_from_fem(df, limit_T, limit_time)

    assert results["max_reached_temperature"].magnitude == 80
    assert results["time_max_reached_temperature"].magnitude == 20
    assert results["constraint_temperature"].magnitude == 0.14285714285714285
    assert results["time_of_demolding"].magnitude == 15


def test_kpi_from_fem_extrapolate1():
    pint.set_application_registry(ureg)  # required to use the same registry
    df = pd.DataFrame(
        {
            "time": pd.Series([0, 10, 20], dtype="pint[hours]"),
            "temperature": pd.Series([10, 40, 80], dtype="pint[degree_Celsius]"),
            "yield": pd.Series([40, 20, 10], dtype="pint[dimensionless]"),
        }
    )

    Q_ = ureg.Quantity
    limit_T = Q_(70, ureg.degC)
    limit_time = 6 * ureg.hours

    results = kpi_from_fem(df, limit_T, limit_time)

    assert results["max_reached_temperature"].magnitude == pytest.approx(130)
    assert results["time_max_reached_temperature"].magnitude == pytest.approx(30)
    assert results["constraint_temperature"].magnitude == pytest.approx(0.8571428571428571)
    assert results["time_of_demolding"].magnitude == pytest.approx(30)


def test_kpi_from_fem_extrapolate2():
    pint.set_application_registry(ureg)  # required to use the same registry
    df = pd.DataFrame(
        {
            "time": pd.Series([0, 10, 20], dtype="pint[hours]"),
            "temperature": pd.Series([10, 40, 80], dtype="pint[degree_Celsius]"),
            "yield": pd.Series([-10, -20, -40], dtype="pint[dimensionless]"),
        }
    )

    Q_ = ureg.Quantity
    limit_T = Q_(70, ureg.degC)
    limit_time = 6 * ureg.hours

    results = kpi_from_fem(df, limit_T, limit_time)

    assert results["max_reached_temperature"].magnitude == pytest.approx(80)
    assert results["time_max_reached_temperature"].magnitude == pytest.approx(20)
    assert results["constraint_temperature"].magnitude == pytest.approx(0.14285714285714285)
    assert results["time_of_demolding"].magnitude == pytest.approx(-10)
