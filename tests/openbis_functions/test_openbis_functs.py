import pytest

from lebedigital.openbis.expstep import ExpStep


def test_load_emodul():
    """Testing if the metadata from yaml emodul files are imported corretly
    """

    yaml_path = '../usecases/MinimumWorkingExample/emodul/metadata_yaml_files/Wolf 8.2 Probe 1.yaml'

    test_sample = ExpStep(name='test_sample')
    test_sample.read_metadata_emodul(yaml_path)
    test_sample.sync_name('name')

    expected_sample_dict = {'name': 'test_sample', 'type': '', 'metadata': {'diameter': 98.6, 'experimentname': 'Wolf 8.2 Probe 1', 'length': 300.2, 'length_unit': 'mm', 'mix_file': '2014_12_10 Wolf.xls', 'operator_date': '28.01.2015', 'operator_timestamp': '10:56:24', 'remark': 'Kraftgeregelt 3,9 kN/s', 'software_specification':
                                                                            'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A', 'specimen_name': 'Wolf 8.2 Probe 1', 'tester_name': 'Kh', 'weight': 5328.0, 'weight_unit': 'g', '$name': 'test_sample'}, 'space': '', 'project': '', 'collection': '', 'parents': [], 'children': [], 'identifier': '', 'permId': '', 'sample_object': None, 'datasets': [], 'dataset_codes': []}

    assert test_sample.__dict__ == expected_sample_dict
