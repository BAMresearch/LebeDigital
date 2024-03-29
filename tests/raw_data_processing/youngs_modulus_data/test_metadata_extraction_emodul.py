from lebedigital.raw_data_processing.youngs_modulus_data.emodul_metadata_extraction import extract_metadata_emodulus
from pathlib import Path

def test_metadata_extraction_emodul():
        """Tesing the metadata extraction on a single example"""

        # setup paths and directories
        data_dir = 'test_data'
        data_path = Path(__file__).parent / data_dir

        mix_file = 'mix.dat'
        specimen_file = 'specimen.dat'

        target_data = {'experimentName': 'test_data',
                                         'software_specification': 'MTS793|MPT|DEU|1|2|,|.|:|49|1|1|A',
                                         'operator_timestamp': '13:25:39',
                                         'operator_date': '01.09.2014',
                                         'tester_name': 'Kh',
                                         'specimen_name': 'BA-Losert E-Modul 28d v. 04.08.14 Probe 4',
                                         'remark': 'Kraftgeregelt 3,9 kN/s',
                                         'weight': 5342.0,
                                         'diameter': 98.6,
                                         'weight_unit': 'g',
                                         'length': 300.3,
                                         'length_unit': 'mm',
                                         'mix_file': '2014_08_05 Rezeptur_MI.xlsx'}

        # run extraction and getting a dictionary with metadata
        test_data = extract_metadata_emodulus(data_path, specimen_file, mix_file)

        assert test_data == target_data
