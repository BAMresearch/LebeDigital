import os
import random
import string
from collections import defaultdict

import pandas as pd
import pytest

from lebedigital.openbis.interbis import Interbis

"""

All tests which require an openBIS login are marked as login

in order to specify which tests to skip these tests run the command

pytest -v -m 'not login' OR pytest -v -m 'login'

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


@pytest.fixture(scope='session', autouse=True)
def o_conn(pytestconfig, sample_code, sample_dict):

    login_val = pytestconfig.getoption('--login')
    password_val = pytestconfig.getoption('--password')

    o = Interbis('https://test.datastore.bam.de/openbis/')

    if login_val != 'no_cl_login' and password_val != 'no_cl_password':
        o.connect_to_datastore(username=login_val, password=password_val)
    else:
        o.connect_to_datastore()

    yield

    created_sample_type = o.get_sample_type(sample_code[0])
    created_sample_type.delete('cleaning up after test run')

    for p_name, p_vals in sample_dict.items():
        prop = o.get_property_type(p_name)
        prop.delete('cleaning up after test run')

    o.logout()


@pytest.mark.login
def test_get_metadata_import_template(o_conn):
    o = Interbis('https://test.datastore.bam.de/openbis/')

    sample_type = 'EXPERIMENTAL_STEP_TEST'

    df = o.get_metadata_import_template(sample_type)

    print(os.getcwd())
    df_expected = pd.read_csv('./gen_import_template.csv', index_col=0, keep_default_na=False)

    pd.testing.assert_frame_equal(df, df_expected)


@pytest.mark.login
def test_get_sample_properties(o_conn):
    o = Interbis('https://test.datastore.bam.de/openbis/')

    sample_type = 'EXPERIMENTAL_STEP_TEST'

    df = o.get_sample_properties(sample_type)

    print(os.getcwd())
    df_expected = pd.read_csv('./gen_sample_properties.csv',
                              index_col=0,
                              keep_default_na=False,
                              )

    df = df.drop(['semanticAnnotations', 'metaData'], axis=1)
    df_expected = df_expected.drop(['semanticAnnotations', 'metaData'], axis=1)

    pd.testing.assert_frame_equal(df, df_expected, check_dtype=False)


@pytest.fixture(scope='session')
def sample_code():
    sample_code = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    sample_prefix = ''.join(sample_code[:5])
    return sample_code, sample_prefix


@pytest.fixture(scope='session')
def sample_dict():
    sample_type_dict = {'TYPESTR': ['VARCHAR', 'typeStr_label', 'typeStr_desc'],
                        'TYPEINT': ['INTEGER', 'typeInt_label', 'typeInt_desc'],
                        'TYPEFLOAT': ['REAL', 'typeFloat_label', 'typeFloat_desc'],
                        'TYPEBOOLEAN': ['BOOLEAN', 'typeBool_label', 'typeBool_desc'],
                        'TYPEHYPERLINK': ['HYPERLINK', 'typeHyperlink_label', 'typeHyperlink_desc'],
                        'TYPEMULTILINEVARCHAR': ['MULTILINE_VARCHAR', 'typeMultilineVarchar_label',
                                                 'typeMultilineVarchar_desc'],
                        'TYPETIMESTAMP': ['TIMESTAMP', 'typeTimestamp_label', 'typeTimestamp_desc']}
    return sample_type_dict


@pytest.mark.login
def test_create_sample_type(sample_code, sample_dict):

    o = Interbis('https://test.datastore.bam.de/openbis/')

    o.create_sample_type(sample_code=sample_code[0], sample_prefix=sample_code[1], sample_properties=sample_dict)

    true_sample_type = o.get_sample_type(sample_code[0]).get_property_assignments()

    true_sample_type_df = true_sample_type.df

    true_sample_props_list = list(true_sample_type_df['propertyType'].values)

    true_sample_dict = defaultdict(list)

    for prop in true_sample_props_list:
        prop_object = o.get_property_type(prop)
        true_sample_dict[prop] = [prop_object.dataType, prop_object.label, prop_object.description]

    assert dict(true_sample_dict) == sample_dict
