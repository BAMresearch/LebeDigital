import pytest

from lebedigital.openbis.expstep import ExpStep
from pybis import Openbis
import pandas as pd

# https://docs.pytest.org/en/7.2.x/how-to/monkeypatch.html

# THE ISSUE HERE IS THAT THE RETURN VALUE FROM THE MOCKED CALL WOULD BE
# A PYBIS THINGS OBJECT AND NOT A DATAFRAME WHICH I NEED TO PARSE STUFF
# FROM IT


class MockResponse:

    # mock get_samples() to return it somehow as the api call
    @staticmethod
    def get_samples():
        pybis_response = {'identifier': {0: '/CKUJATH/TEST_AMBETON/EXP_VISKO867278', 1: '/CKUJATH/TEST_AMBETON/EXP_VISKO867295', 2: '/CKUJATH/TEST_AMBETON/EXP_VISKO867292', 3: '/CKUJATH/TEST_AMBETON/EXP_VISKO867291', 4: '/CKUJATH/TEST_AMBETON/EXP_VISKO867294', 5: '/CKUJATH/TEST_AMBETON/EXP_VISKO867293', 6: '/CKUJATH/TEST_AMBETON/EXP_VISKO867288', 7: '/CKUJATH/TEST_AMBETON/EXP_VISKO867287', 8: '/CKUJATH/TEST_AMBETON/EXP_VISKO867290', 9: '/CKUJATH/TEST_AMBETON/EXP_VISKO867289', 10: '/CKUJATH/TEST_AMBETON/EXP_VISKO867284', 11: '/CKUJATH/TEST_AMBETON/EXP_VISKO867283', 12: '/CKUJATH/TEST_AMBETON/EXP_VISKO867286', 13: '/CKUJATH/TEST_AMBETON/EXP_VISKO867285', 14: '/CKUJATH/TEST_AMBETON/EXP_VISKO867280', 15: '/CKUJATH/TEST_AMBETON/EXP_VISKO867279', 16: '/CKUJATH/TEST_AMBETON/EXP_VISKO867282', 17: '/CKUJATH/TEST_AMBETON/EXP_VISKO867281'}, 'permId': {0: '20221018105358763-867795', 1: '20221018105502129-867829', 2: '20221018105451295-867823', 3: '20221018105447690-867821', 4: '20221018105456578-867827', 5: '20221018105454270-867825', 6: '20221018105437682-867815', 7: '20221018105434201-867813', 8: '20221018105444081-867819', 9: '20221018105440934-867817', 10: '20221018105420025-867807', 11: '20221018105417357-867805', 12: '20221018105430194-867811', 13: '20221018105424950-867809', 14: '20221018105406726-867799', 15: '20221018105403230-867797', 16: '20221018105413920-867803', 17: '20221018105410053-867801'}, 'type': {0: 'EXPERIMENTAL_STEP_VISKO', 1: 'EXPERIMENTAL_STEP_VISKO', 2: 'EXPERIMENTAL_STEP_VISKO', 3: 'EXPERIMENTAL_STEP_VISKO', 4: 'EXPERIMENTAL_STEP_VISKO', 5: 'EXPERIMENTAL_STEP_VISKO', 6: 'EXPERIMENTAL_STEP_VISKO', 7: 'EXPERIMENTAL_STEP_VISKO', 8: 'EXPERIMENTAL_STEP_VISKO', 9: 'EXPERIMENTAL_STEP_VISKO', 10: 'EXPERIMENTAL_STEP_VISKO', 11: 'EXPERIMENTAL_STEP_VISKO', 12: 'EXPERIMENTAL_STEP_VISKO', 13: 'EXPERIMENTAL_STEP_VISKO', 14: 'EXPERIMENTAL_STEP_VISKO', 15: 'EXPERIMENTAL_STEP_VISKO', 16: 'EXPERIMENTAL_STEP_VISKO', 17: 'EXPERIMENTAL_STEP_VISKO'}, 'registrator': {0: 'araderma', 1: 'araderma', 2: 'araderma', 3: 'araderma', 4: 'araderma', 5: 'araderma', 6: 'araderma', 7: 'araderma',
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          8: 'araderma', 9: 'araderma', 10: 'araderma', 11: 'araderma', 12: 'araderma', 13: 'araderma', 14: 'araderma', 15: 'araderma', 16: 'araderma', 17: 'araderma'}, 'registrationDate': {0: '2022-10-18 10:53:59', 1: '2022-10-18 10:55:02', 2: '2022-10-18 10:54:51', 3: '2022-10-18 10:54:48', 4: '2022-10-18 10:54:57', 5: '2022-10-18 10:54:54', 6: '2022-10-18 10:54:38', 7: '2022-10-18 10:54:34', 8: '2022-10-18 10:54:44', 9: '2022-10-18 10:54:41', 10: '2022-10-18 10:54:20', 11: '2022-10-18 10:54:17', 12: '2022-10-18 10:54:30', 13: '2022-10-18 10:54:25', 14: '2022-10-18 10:54:07', 15: '2022-10-18 10:54:03', 16: '2022-10-18 10:54:14', 17: '2022-10-18 10:54:10'}, 'modifier': {0: 'araderma', 1: 'araderma', 2: 'araderma', 3: 'araderma', 4: 'araderma', 5: 'araderma', 6: 'araderma', 7: 'araderma', 8: 'araderma', 9: 'araderma', 10: 'araderma', 11: 'araderma', 12: 'araderma', 13: 'araderma', 14: 'araderma', 15: 'araderma', 16: 'araderma', 17: 'araderma'}, 'modificationDate': {0: '2022-10-18 11:49:11', 1: '2022-10-18 13:23:10', 2: '2022-10-18 13:05:19', 3: '2022-10-18 13:05:26', 4: '2022-10-18 13:06:17', 5: '2022-10-18 13:06:21', 6: '2022-10-18 13:04:16', 7: '2022-10-18 13:04:22', 8: '2022-10-18 13:05:33', 9: '2022-10-18 13:04:27', 10: '2022-10-18 13:02:03', 11: '2022-10-18 13:00:04', 12: '2022-10-18 13:02:10', 13: '2022-10-18 13:02:19', 14: '2022-10-18 11:51:08', 15: '2022-10-18 11:51:19', 16: '2022-10-18 13:00:12', 17: '2022-10-18 13:00:18'}, '$NAME': {0: '3Dm3_0_1rpm_Vogel_2_7_T17_01', 1: '3Dm3_0_1rpm_Vogel_2_7_T42_03', 2: '3Dm3_0_1rpm_Vogel_2_7_T37_03', 3: '3Dm3_0_1rpm_Vogel_2_7_T37_02', 4: '3Dm3_0_1rpm_Vogel_2_7_T42_02', 5: '3Dm3_0_1rpm_Vogel_2_7_T42_01', 6: '3Dm3_0_1rpm_Vogel_2_7_T32_02', 7: '3Dm3_0_1rpm_Vogel_2_7_T32_01', 8: '3Dm3_0_1rpm_Vogel_2_7_T37_01', 9: '3Dm3_0_1rpm_Vogel_2_7_T32_03', 10: '3Dm3_0_1rpm_Vogel_2_7_T27_01', 11: '3Dm3_0_1rpm_Vogel_2_7_T22_03', 12: '3Dm3_0_1rpm_Vogel_2_7_T27_03', 13: '3Dm3_0_1rpm_Vogel_2_7_T27_02', 14: '3Dm3_0_1rpm_Vogel_2_7_T17_03', 15: '3Dm3_0_1rpm_Vogel_2_7_T17_02', 16: '3Dm3_0_1rpm_Vogel_2_7_T22_02', 17: '3Dm3_0_1rpm_Vogel_2_7_T22_01'}}
        mock_df = pd.DataFrame.from_dict(pybis_response)
        return mock_df


def test_get_sample_names(monkeypatch):

    # Any arguments may be passed and mock_get() will always return our
    # mocked object, which only has the .json() method.
    def mock_get(*args, **kwargs):
        return MockResponse()

    # apply the monkeypatch for requests.get to mock_get
    monkeypatch.setattr(Openbis, "get_samples", mock_get)

    o = Openbis("https://fakeurl")
    space = 'CKUJATH'
    project = 'TEST_AMBETON'
    collection = '/CKUJATH/TEST_AMBETON/VISKO_DATA_COLLECTION'

    # app.get_json, which contains requests.get, uses the monkeypatch
    result = ExpStep.get_sample_names(o, space, project, collection)

    expected_result = ['EXP_VISKO867278 (3Dm3_0_1rpm_Vogel_2_7_T17_01)', 'EXP_VISKO867295 (3Dm3_0_1rpm_Vogel_2_7_T42_03)', 'EXP_VISKO867292 (3Dm3_0_1rpm_Vogel_2_7_T37_03)', 'EXP_VISKO867291 (3Dm3_0_1rpm_Vogel_2_7_T37_02)', 'EXP_VISKO867294 (3Dm3_0_1rpm_Vogel_2_7_T42_02)', 'EXP_VISKO867293 (3Dm3_0_1rpm_Vogel_2_7_T42_01)', 'EXP_VISKO867288 (3Dm3_0_1rpm_Vogel_2_7_T32_02)', 'EXP_VISKO867287 (3Dm3_0_1rpm_Vogel_2_7_T32_01)', 'EXP_VISKO867290 (3Dm3_0_1rpm_Vogel_2_7_T37_01)',
                       'EXP_VISKO867289 (3Dm3_0_1rpm_Vogel_2_7_T32_03)', 'EXP_VISKO867284 (3Dm3_0_1rpm_Vogel_2_7_T27_01)', 'EXP_VISKO867283 (3Dm3_0_1rpm_Vogel_2_7_T22_03)', 'EXP_VISKO867286 (3Dm3_0_1rpm_Vogel_2_7_T27_03)', 'EXP_VISKO867285 (3Dm3_0_1rpm_Vogel_2_7_T27_02)', 'EXP_VISKO867280 (3Dm3_0_1rpm_Vogel_2_7_T17_03)', 'EXP_VISKO867279 (3Dm3_0_1rpm_Vogel_2_7_T17_02)', 'EXP_VISKO867282 (3Dm3_0_1rpm_Vogel_2_7_T22_02)', 'EXP_VISKO867281 (3Dm3_0_1rpm_Vogel_2_7_T22_01)']

    assert result == expected_result
