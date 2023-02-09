import os
import random
import string
from datetime import datetime, timedelta
from enum import Enum

import pandas as pd
import pytest

from lebedigital.openbis.interbis import Interbis

"""

All tests which require an openBIS login are marked as login

in order to specify which tests to skip these tests run the command

pytest -m 'not login' OR pytest -m 'login' --login <username> --password <password> in test folder

"""


class Constants(Enum):
    space: str = 'DUMMY'
    project: str = 'TEST_PROJECT'
    collection: str = 'TEST_COLLECTION'
    collection_id: str = '/DUMMY/TEST_PROJECT/TEST_COLLECTION'
    sample_type: str = 'EXPERIMENTAL_STEP'
    db_url: str = "https://openbis.matolab.org/openbis/"
    testing_sample_name: str = 'TESTING_SAMPLE_NAME_PYTEST_DO_NOT_DELETE'


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


@pytest.fixture(scope='session', autouse=True)
def setup(pytestconfig):
    login_val = pytestconfig.getoption('--login')
    password_val = pytestconfig.getoption('--password')

    o = Interbis(Constants.db_url.value)

    if login_val != 'no_cl_login' and password_val != 'no_cl_password':
        o.connect_to_datastore(username=login_val, password=password_val)
    else:
        o.connect_to_datastore()

    # Create project and collection, if not there
    try:
        o.get_project(projectId=f"/{Constants.space.value}/{Constants.project.value}")
    except ValueError:
        project_obj = o.new_project(space=Constants.space.value, code=Constants.project.value, description="Test project")
        project_obj.save()

    try:
        o.get_collection(code=Constants.collection_id.value)
    except ValueError:
        collection_obj = o.new_collection(project=Constants.project.value, code=Constants.collection.value, type="COLLECTION")
        collection_obj.save()

    # Creating testing object
    sample = o.new_sample(
        type=Constants.sample_type.value,
        space=Constants.space.value,
        project=Constants.project.value,
        collection=Constants.collection_id.value
    )

    sample.set_props({
        '$name': Constants.testing_sample_name.value,
        'start_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'end_date': (datetime.now() + timedelta(hours=2)).strftime("%Y-%m-%d %H:%M:%S"),
        'experimental_step.experimental_goals': 'Testing',
        'experimental_step.experimental_description': 'Also testing',
    })

    sample.save()

    yield

    sample.delete('cleaning up after test run')

    o.logout()


@pytest.mark.login
def test_get_metadata_import_template(setup):
    o = Interbis(Constants.db_url.value)

    df = o.get_metadata_import_template(Constants.sample_type.value)

    df_expected = pd.read_csv('./gen_import_template.csv', index_col=0, keep_default_na=False)

    pd.testing.assert_frame_equal(df, df_expected)


@pytest.mark.login
def test_get_sample_type_properties(setup):
    o = Interbis(Constants.db_url.value)

    df = o.get_sample_type_properties(Constants.sample_type.value)

    df_expected = pd.read_csv('./gen_sample_properties.csv',
                              index_col=0,
                              keep_default_na=False,
                              )

    df = df.drop(['semanticAnnotations', 'metaData'], axis=1)
    df_expected = df_expected.drop(['semanticAnnotations', 'metaData'], axis=1)

    pd.testing.assert_frame_equal(df, df_expected, check_dtype=False)


# no way of testing on dummy account
@pytest.mark.skip
def test_create_sample_type(sample_code, sample_dict):
    o = Interbis(Constants.db_url.value)

    o.create_sample_type(sample_code=sample_code[0], sample_prefix=sample_code[1], sample_properties=sample_dict)

    true_sample_type = o.get_sample_type(sample_code[0]).get_property_assignments()

    true_sample_type_df = true_sample_type.df

    true_sample_props_list = list(true_sample_type_df['propertyType'].values)

    true_sample_dict = {}

    for prop in true_sample_props_list:
        prop_object = o.get_property_type(prop)
        true_sample_dict[prop] = [prop_object.dataType, prop_object.label, prop_object.description]

    assert true_sample_dict == sample_dict


@pytest.mark.login
@pytest.mark.parametrize("sample, output", [(Constants.testing_sample_name.value, True),
                                            (''.join(random.choice('0123456789ABCDEF') for _ in range(16)), False)])
def test_exists_in_datastore(setup, sample, output):
    o = Interbis(Constants.db_url.value)

    assert o.exists_in_datastore(str(sample)) == output
