import os
from lebedigital.raw_data_processing.mixture \
    .mixture_metadata_extraction import extract_metadata_mixture
from config import Config
import unittest


class TestMixtureMetadataExtraction(unittest.TestCase):
    def test_extraction(self):
        """Testing the mixture metadata extraction on a single example"""

        # setting up the test example input
        input_path = '../usecases/MinimumWorkingExample/Data/Mischungen/2014_08_05 Rezeptur_MI.xlsx'

        target_data = Config.TEST_MIXTURE_METADATA_EXTRACTION_EXPECTED_VALUE
        # run extraction and getting a dictionary with metadata
        test_data = extract_metadata_mixture(input_path, None)

        # make sure if the path exists
        self.assertTrue(os.path.exists(input_path))

        # checking if result is correct
        self.assertEqual(test_data.__str__(), target_data)


if __name__ == '__main__':
    unittest.main()
