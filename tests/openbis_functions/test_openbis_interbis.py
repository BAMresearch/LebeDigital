import pandas as pd
import pytest
import random, string

from lebedigital.openbis.interbis import Interbis

"""

All tests which require an openBIS login are marked as login

in order to specify which tests to skip these tests run the command

pytest -v -m

"""


def connect_for_tests(login: str, password: str) -> Interbis:
    """ Different login for tests, allows to parse login/password as both
    command line arguments and while the tests are running

    Args:
        login (str): Login for openBIS
        password (str): Password for openBIS

    Returns:
        Interbis: connected Interbis object
    """

    o = Interbis('https://test.datastore.bam.de/openbis/')
    if login != 'no_cl_login' and password != 'no_cl_password':
        o.connect_to_datastore(username=login, password=password)
    else:
        o.connect_to_datastore()

    return o


@pytest.mark.login
def test_get_metadata_import_template(login, password):
    o = connect_for_tests(login, password)

    sample_type = 'EXPERIMENTAL_STEP_TEST'

    df = o.get_metadata_import_template(sample_type)

    df_expected = pd.read_csv('./gen_import_template.csv', index_col=0, keep_default_na=False)

    o.logout()

    pd.testing.assert_frame_equal(df, df_expected)


@pytest.mark.login
def test_get_sample_properties(login, password):
    o = connect_for_tests(login, password)

    sample_type = 'EXPERIMENTAL_STEP_TEST'

    df = o.get_sample_properties(sample_type)

    df_expected = pd.read_csv('./gen_sample_properties.csv',
                              index_col=0,
                              keep_default_na=False,
                              )

    df = df.drop(['semanticAnnotations', 'metaData'], axis=1)
    df_expected = df_expected.drop(['semanticAnnotations', 'metaData'], axis=1)

    o.logout()

    pd.testing.assert_frame_equal(df, df_expected, check_dtype=False)


@pytest.mark.login
def test_create_sample_type(login, password):
    o = connect_for_tests(login, password)

    sample_code = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    sample_prefix = ''.join(sample_code[:5])

    sample_type_dict = {
        'typeStr': ['VARCHAR', 'typeStr_label', 'typeStr_desc'],
        'typeInt': ['INTEGER', 'typeInt_label', 'typeInt_desc'],
        'typeFloat': ['REAL', 'typeFloat_label', 'typeFloat_desc'],
        'typeBoolean': ['BOOLEAN', 'typeBool_label', 'typeBool_desc'],
    }

    o.create_sample_type()

    o.logout()
