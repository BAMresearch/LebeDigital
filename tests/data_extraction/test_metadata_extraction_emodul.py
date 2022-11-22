import unittest
from config import Config

from lebedigital.raw_data_processing.youngs_modulus_data \
    .emodul_metadata_extraction import extract_metadata_emodulus
class TestMetaDataExtractionEmodul(unittest.TestCase):
    def test_metadata_extraction_emodul(self):
        """Tesing the metadata extraction on a single example"""

        # setting up the test example
        input = '../usecases/MinimumWorkingExample/Data/E-modul/BA-Losert MI E-Modul 28d v. 04.08.14 Probe 4'
        mix_file = 'mix.dat'
        specimen_file = 'specimen.dat'

        target_data = Config.TEST_METADATA_EXTRACTION_EXPECTED_VALUE

        # run extraction and getting a dictionary with metadata
        test_data = extract_metadata_emodulus(input, specimen_file, mix_file)

        # checking if result is correct
        self.assertEqual(test_data, target_data)