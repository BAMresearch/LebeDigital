import random
import string
from datetime import datetime, timedelta
from enum import Enum
import json
from pathlib import Path
import os

import pandas as pd
import pytest

from lebedigital.openbis.interbis import Interbis

"""

All tests which require an openBIS login are marked as login

in order to specify which tests to skip these tests run the command

pytest -m 'not login' OR pytest -m 'login' --login <username> --password <password> --url <url> in test folder

url is localhost by default, when you want to test on some other database then you need to specify the url address
"""


class Constants(Enum):
    space: str = 'DEFAULT'
    project: str = 'TEST_PROJECT'
    collection: str = 'TEST_COLLECTION'
    collection_id: str = '/DEFAULT/TEST_PROJECT/TEST_COLLECTION'
    sample_type: str = 'EXPERIMENTAL_STEP'
    # db_url: str = "https://openbis.matolab.org/openbis/"
    db_url: str = "https://localhost:8443/openbis/"
    testing_sample_name: str = 'PYTEST_SAMPLE'
    testing_sample_identifier: str = "/DEFAULT/TEST_PROJECT/PYTEST_SAMPLE"
    parent_hint_label: str = "testing label"
    sample_type_typechecker_code: str = "ST_TYPECHECKER"
    sample_type_typechecker_prefix: str = "TYPECHK"


class Filepaths(Enum):
    import_template: Path = Path('./openbis_functions/test_files/gen_import_template.csv')
    sample_properties: Path = Path('./openbis_functions/test_files/gen_sample_properties.csv')
    test_sheet: Path = Path('./openbis_functions/test_files/test_sheet.xlsx')
    init_settings: Path = Path('./openbis_functions/test_files/init_settings.json')
    excel_output: Path = Path('./openbis_functions/test_files/output.xlsx')
    filled_out_sheet: Path = Path('./openbis_functions/test_files/filled_out_sheet.xlsx')


test_results = {
    "idx_parent_hint": None
}

created_samples_in_tests = []


@pytest.fixture(scope='session')
def sample_code():
    sample_code = ''.join(random.choices(
        string.ascii_letters + string.digits, k=10))
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


@pytest.fixture(scope='session')
def expected_df_import():
    return pd.read_csv(Filepaths.import_template.value, index_col=0, keep_default_na=False)


@pytest.fixture(scope='session', autouse=True)
def setup(pytestconfig):
    login_val = pytestconfig.getoption('--login')
    password_val = pytestconfig.getoption('--password')
    chosen_runner = pytestconfig.getoption('--url')

    # o = Interbis(Constants.db_url.value)
    o = Interbis(chosen_runner, verify_certificates=False)

    if login_val != 'no_cl_login' and password_val != 'no_cl_password':
        o.connect_to_datastore(username=login_val, password=password_val)
    else:
        o.connect_to_datastore()

    # Setting initial settings to make parent_hint setting possible
    with open(Filepaths.init_settings.value, 'r') as file:
        default_settings = json.load(file)

    settings_sample = o.get_sample("/ELN_SETTINGS/GENERAL_ELN_SETTINGS")
    settings_sample.props["$eln_settings"] = json.dumps(default_settings)
    settings_sample.save()

    # Create project and collection, if not there
    try:
        project_obj = o.get_project(
            projectId=f"/{Constants.space.value}/{Constants.project.value}")
    except (KeyError, ValueError):
        project_obj = o.new_project(space=Constants.space.value, code=Constants.project.value,
                                    description="Test project")
        project_obj.save()

    try:
        collection_obj = o.get_collection(code=Constants.collection_id.value)
    except ValueError:
        collection_obj = o.new_collection(project=Constants.project.value, code=Constants.collection.value,
                                          type="COLLECTION")
        collection_obj.save()

    # Creating testing object
    sample = o.new_sample(
        type=Constants.sample_type.value,
        space=Constants.space.value,
        project=Constants.project.value,
        collection=Constants.collection_id.value,
        code=Constants.testing_sample_name.value,
    )

    sample.set_props({
        '$name': Constants.testing_sample_name.value,
        'start_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'end_date': (datetime.now() + timedelta(hours=2)).strftime("%Y-%m-%d %H:%M:%S"),
        'experimental_step.experimental_goals': 'Testing',
        'experimental_step.experimental_description': 'Also testing',
    })

    sample.save()
    created_samples_in_tests.append(sample.identifier)

    test_vocab_terms = [
        {"code": 'very', "label": "very", "description": "very"},
        {"code": 'interesting', "label": "interesting", "description": "interesting"},
        {"code": 'terms', "label": "terms", "description": "terms"},
        {"code": 1, "label": "1", "description": "1"},
        {"code": 0.25, "label": "0.25", "description": "0.25"},
    ]

    test_vocab = o.new_vocabulary(
        code='test_vocab',
        description='test_vocab_description',
        terms=test_vocab_terms
    )

    test_vocab.save()

    typechecker_sample_type_props = {
        'testing_vocabulary': [
            'CONTROLLEDVOCABULARY',
            'testing_vocabulary_label',
            'testing_vocabulary_description',
            'test_vocab'
        ],
        'testing_real': [
            'REAL',
            'testing_real_label',
            'testing_real_description',
        ],
        'testing_timestamp': [
            'TIMESTAMP',
            'testing_timestamp_label',
            'testing_timestamp_description',
        ],
        'testing_varchar': [
            'VARCHAR',
            'testing_varchar_label',
            'testing_varchar_description',
        ]
    }

    o.create_sample_type(
        Constants.sample_type_typechecker_code.value,
        Constants.sample_type_typechecker_prefix.value,
        typechecker_sample_type_props
    )

    yield

    for sample_id in created_samples_in_tests:
        sample = o.get_sample(sample_id)
        sample.delete("cleaning up after tests")

    o.logout()


@pytest.mark.login
@pytest.mark.parametrize("write, sheet_name, path", [(False, "metadata", ""),
                                                     (True, "some_named_sheet", Filepaths.excel_output.value)])
def test_get_metadata_import_template(setup, pytestconfig, expected_df_import, write, sheet_name, path):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    df = o.get_metadata_import_template(Constants.sample_type.value, write, sheet_name, path)

    if write:
        os.path.isfile(path)
        written_df = pd.read_excel(path, sheet_name=sheet_name, index_col=0, keep_default_na=False)
        pd.testing.assert_frame_equal(written_df, expected_df_import)
        os.remove(path)
    else:
        pd.testing.assert_frame_equal(df, expected_df_import)


@pytest.mark.login
def test_import_props_from_template(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    read_sample_props = o.import_props_from_template(Filepaths.filled_out_sheet.value)

    checked_keys = ["experimental_step.experimental_goals", "experimental_step.experimental_description", "$name"]

    expected_sample_props = o.get_sample(Constants.testing_sample_identifier.value).props()
    expected_sample_props = {
        key: val
        for key, val in expected_sample_props.items()
        if key in checked_keys
    }

    assert read_sample_props == expected_sample_props


@pytest.mark.login
def test_get_sample_dict(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    sample_dict = o.get_sample_dict(Constants.testing_sample_identifier.value)

    # Dict with:
    #   mentions of date removed
    #   removed permId

    expected_sample_dict = {
        '$ANNOTATIONS_STATE': None,
        '$NAME': 'PYTEST_SAMPLE',
        '$SHOW_IN_PROJECT_OVERVIEW': None,
        '$XMLCOMMENTS': None,
        'EXPERIMENTAL_STEP.EXPERIMENTAL_DESCRIPTION': 'Also testing',
        'EXPERIMENTAL_STEP.EXPERIMENTAL_GOALS': 'Testing',
        'EXPERIMENTAL_STEP.EXPERIMENTAL_RESULTS': None,
        'EXPERIMENTAL_STEP.SPREADSHEET': None,
        'FINISHED_FLAG': None,
        'NOTES': None,
        'PUBLICATION': None,
        'REFERENCE': None,
        'children': [],
        'code': 'PYTEST_SAMPLE',
        'collection': '/DEFAULT/TEST_PROJECT/TEST_COLLECTION',
        'components': [],
        'container': None,
        'experiment': '/DEFAULT/TEST_PROJECT/TEST_COLLECTION',
        'frozen': False,
        'frozenForChildren': False,
        'frozenForComponents': False,
        'frozenForDataSets': False,
        'frozenForParents': False,
        'identifier': '/DEFAULT/TEST_PROJECT/PYTEST_SAMPLE',
        'modifier': 'admin',
        'parents': [],
        'project': '/DEFAULT/TEST_PROJECT',
        'registrator': 'admin',
        'space': 'DEFAULT',
        'tags': [],
        'type': 'EXPERIMENTAL_STEP'
    }

    assert expected_sample_dict.items() <= sample_dict.items()


@pytest.mark.login
@pytest.mark.parametrize('level', [('full'), ('space'), ('project'), ('collection')])
def test_get_overview(setup, pytestconfig, level):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    overview = o.get_overview(level=level, space=Constants.space.value, project=Constants.project.value, collection=Constants.collection.value)

    print('\n======================== ' + level.upper() + ' ========================')
    print(overview)

    if level == 'full':
        assert overview['DATASTORE']['ELN_SETTINGS']
    elif level == 'space':
        assert overview['DEFAULT']['TEST_PROJECT']
    elif level == 'project':
        assert overview['TEST_PROJECT']['TEST_COLLECTION']
    else:
        assert len(overview['TEST_COLLECTION']) >= 1


@pytest.mark.login
def test_get_sample_type_properties(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    df = o.get_sample_type_properties(Constants.sample_type.value)

    df_expected = pd.read_csv(Filepaths.sample_properties.value,
                              index_col=0,
                              keep_default_na=False,
                              )

    ignored_columns = ['semanticAnnotations', 'metaData', 'registrationDate']

    df = df.drop(ignored_columns, axis=1)
    df_expected = df_expected.drop(ignored_columns, axis=1)

    pd.testing.assert_frame_equal(df, df_expected, check_dtype=False)


@pytest.mark.login
def test_create_sample_type(sample_code, sample_dict, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    o.create_sample_type(
        sample_code=sample_code[0], sample_prefix=sample_code[1], sample_properties=sample_dict)

    true_sample_type = o.get_sample_type(
        sample_code[0]).get_property_assignments()

    true_sample_type_df = true_sample_type.df

    true_sample_props_list = list(true_sample_type_df['propertyType'].values)

    true_sample_dict = {}

    for prop in true_sample_props_list:
        prop_object = o.get_property_type(prop)
        true_sample_dict[prop] = [prop_object.dataType,
                                  prop_object.label, prop_object.description]

    assert true_sample_dict == sample_dict

    o.get_sample_type(sample_code[0]).delete(reason='testing cleanup')


@pytest.mark.login
@pytest.mark.parametrize("sample_name, output", [(Constants.testing_sample_name.value, True),
                                                 (''.join(random.choice('0123456789ABCDEF') for _ in range(16)), False),
                                                 ("make_two_of_me", True)])
def test_exists_in_datastore(setup, sample_name, output, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    if sample_name == "make_two_of_me":
        for _ in range(2):
            new_sample = o.new_sample(
                type=Constants.sample_type.value,
                space=Constants.space.value,
                project=Constants.project.value,
                collection=Constants.collection_id.value
            )

            new_sample.set_props({
                '$name': sample_name,
            })

            new_sample.save()
            created_samples_in_tests.append(new_sample.identifier)

    assert o.exists_in_datastore(str(sample_name)) == output


@pytest.mark.login
def test_get_sample_identifier(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    fetched_sample_identifier = o.get_sample_identifier(Constants.testing_sample_name.value)

    assert fetched_sample_identifier == f"/{Constants.space.value}/{Constants.project.value}/{Constants.testing_sample_name.value}"


@pytest.mark.login
@pytest.mark.parametrize("collection, should_pass",
                         [(Constants.collection.value, True),
                          (''.join(random.choice('0123456789ABCDEF') for _ in range(16)), False)])
def test_get_collection_identifier(setup, pytestconfig, collection, should_pass):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    if should_pass:
        collection_identifier = o.get_collection_identifier(collection)
        assert collection_identifier == Constants.collection_id.value
    else:
        with pytest.raises(ValueError) as e_info:
            collection_identifier = o.get_collection_identifier(collection)


@pytest.mark.login
def test_create_parent_hint(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    o.create_parent_hint(sample_type=Constants.sample_type.value, label="testing label", parent_type=Constants.sample_type.value)

    settings_sample = o.get_sample("/ELN_SETTINGS/GENERAL_ELN_SETTINGS")

    settings = json.loads(settings_sample.props["$eln_settings"])

    output = settings["sampleTypeDefinitionsExtension"][Constants.sample_type.value]["SAMPLE_PARENTS_HINT"]
    output = [idx for idx, val in enumerate(output) if val["LABEL"] == Constants.parent_hint_label.value]

    test_results["idx_parent_hint"] = output

    assert len(output) == 1


@pytest.mark.login
def test_set_parent_hint(setup, pytestconfig):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    comment_value = 'comment_comment'

    parent_sample = o.new_sample(
        type=Constants.sample_type.value,
        space=Constants.space.value,
        project=Constants.project.value,
        collection=Constants.collection_id.value,
        children=[Constants.testing_sample_identifier.value]
    )

    parent_sample.set_props({
        '$name': "parent_sample_hint_test",
    })

    parent_sample.save()
    created_samples_in_tests.append(parent_sample.identifier)

    o.set_parent_annotation(
        child_sample=Constants.testing_sample_identifier.value,
        parent_sample=parent_sample.identifier,
        comment=comment_value,
    )

    annotation = o.get_parent_annotations(Constants.testing_sample_identifier.value)

    assert annotation[parent_sample.permId]['parentAnnotations']['ANNOTATION.SYSTEM.COMMENTS'] == comment_value


@pytest.mark.login
@pytest.mark.parametrize("param_name, param_val",
                         [('testing_vocabulary', "interesting"),
                          ('testing_vocabulary', 1),
                          ('testing_vocabulary', 0.25)])
def test_generate_typechecker_passing(setup, pytestconfig, param_name, param_val):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    sample = o.new_sample(
        type=Constants.sample_type_typechecker_code.value,
        space=Constants.space.value,
        project=Constants.project.value,
        collection=Constants.collection_id.value,
    )

    sample_props = {
        'testing_varchar': 'varchar',
        'testing_real': '21.37',
        'testing_timestamp': '10.05.2023 10:05',
    }

    sample_props = sample_props | {param_name: param_val}

    Model = o.generate_typechecker(Constants.sample_type_typechecker_code.value)

    model_return = Model(**sample_props)

    assert model_return.testing_varchar == 'varchar'
    assert model_return.testing_real == 21.37
    assert model_return.testing_timestamp == '2023-10-05 10:05'
    assert model_return.testing_vocabulary == str(param_val).upper()

    sample_props = model_return.dict(exclude_unset=True)

    sample.set_props(sample_props)

    sample.save()

    created_samples_in_tests.append(sample.identifier)


@ pytest.mark.login
@ pytest.mark.xfail
@ pytest.mark.parametrize("param_name, param_val",
                          [('testing_timestamp', 'not_a_date'),
                           ('testing_vocabulary', '🤨'),
                           ('testing_real', 'cant_cast_this')])
def test_generate_typechecker_failing(setup, pytestconfig, param_name, param_val):

    chosen_runner = pytestconfig.getoption('--url')
    o = Interbis(chosen_runner, verify_certificates=False)

    sample_props = {param_name: param_val}

    Model = o.generate_typechecker(Constants.sample_type_typechecker_code.value)

    Model(**sample_props)
    # should fail here
    # should fail here
    # should fail here
