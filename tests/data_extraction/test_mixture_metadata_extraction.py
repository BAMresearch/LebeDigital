import os
from lebedigital.raw_data_processing.mixture \
    .mixture_metadata_extraction import extract_metadata_mixture
import unittest


class TestMixtureMetadataExtraction(unittest.TestCase):
    def test_extraction(self):
        """Testing the mixture metadata extraction on a single example"""

        # setting up the test example input
        input_path = '../usecases/MinimumWorkingExample/Data/Mischungen/2014_08_05 Rezeptur_MI.xlsx'

        target_data = "{'operator_date': '2014-06-30', 'tester_name': 'Werner', 'specimen_name': 'BA-Losert MI', 'cement--QuantityInMix': 330.0, 'cement--BulkDensity': 3.123, 'cement--Volume': 105.7, 'cement--Annotation': 'CEM I 42.5 R', 'water_total--QuantityInMix': 175.0, 'water_total--BulkDensity': 1.0, 'water_total--Volume': 175.0, 'water_cement_ratio': 0.5303030303030303, 'water_effective--QuantityInMix': nan, 'water_effective--BulkDensity': nan, 'water_effective--Volume': nan, 'air_content--QuantityInMix': 0.0, 'air_content--BulkDensity': 0.0, 'air_content--Volume': 20.0, 'addition1--QuantityInMix': 273.0, 'addition1--BulkDensity': 2.74, 'addition1--Volume': 99.6, 'addition1--Annotation': 'Medenbach - Kalksteinmehl', 'admixture--QuantityInMix': 5.61, 'admixture--BulkDensity': 1.05, 'admixture--Volume': 5.3, 'admixture--Annotation': 'FM 595 BASF', 'aggregate--QuantityInMix': 1564.0, 'aggregate--BulkDensity': nan, 'aggregate--Volume': 594.4000000000001}"
        # run extraction and getting a dictionary with metadata
        test_data = extract_metadata_mixture(input_path, None)

        # make sure if the path exists
        self.assertTrue(os.path.exists(input_path))

        # checking if result is correct
        self.assertEqual(test_data.__str__(), target_data)


if __name__ == '__main__':
    unittest.main()
